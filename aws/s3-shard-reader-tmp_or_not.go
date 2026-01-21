package main

import (
	"context"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"
	"github.com/aws/aws-sdk-go-v2/feature/s3/manager"
	"github.com/aws/aws-sdk-go-v2/service/s3"
)

const (
	DefaultChunkSize = 8 * 1024 * 1024 // 8MB for S3 parts
)

func main() {
	// --- 1. Argument Parsing ---
	if len(os.Args) < 4 {
		log.Fatalf("Usage: %s <s3_key> <output_fifo> <byte_range> [key=val...]", os.Args[0])
	}

	s3Key := os.Args[1]
	outputFifoPath := os.Args[2]
	byteRangeStr := os.Args[3]

	kvArgs := parseArgs(os.Args[4:])
	bucket := os.Getenv("AWS_BUCKET")
	if bucket == "" {
		log.Fatal("Error: AWS_BUCKET environment variable not set")
	}

	debug := kvArgs["debug"] == "true"
	
	// Determine RAM buffer threshold (default 128MB)
	ramLimitMB := 128
	if envLimit := os.Getenv("PASH_S3_RAM_BUFFER_MB"); envLimit != "" {
		if val, err := strconv.Atoi(envLimit); err == nil {
			ramLimitMB = val
		}
	}
	ramLimitBytes := int64(ramLimitMB) * 1024 * 1024

	// --- 2. Setup S3 ---
	startByte, endByte, err := parseByteRange(byteRangeStr)
	if err != nil {
		log.Fatalf("Invalid byte range: %v", err)
	}
	contentLength := endByte - startByte + 1

	if debug {
		log.Printf("[GO-MAIN] S3: %s/%s | Size: %.2f MB", bucket, s3Key, float64(contentLength)/1024/1024)
	}

	ctx := context.TODO()
	cfg, err := config.LoadDefaultConfig(ctx)
	if err != nil {
		log.Fatalf("unable to load SDK config, %v", err)
	}
	client := s3.NewFromConfig(cfg)

	// --- 3. Execution Strategy ---
	
	// Strategy A: IN-MEMORY (Fastest)
	// If shard size fits in allowed RAM, skip disk entirely.
	if contentLength <= ramLimitBytes {
		if debug {
			log.Printf("[GO-SMART] Strategy: RAM BUFFER (Size < %dMB)", ramLimitMB)
		}
		processInMemory(ctx, client, bucket, s3Key, startByte, endByte, outputFifoPath, debug)
	} else {
		// Strategy B: DISK BUFFER (Safe)
		// If shard is huge, use temp file to avoid OOM.
		if debug {
			log.Printf("[GO-SMART] Strategy: DISK BUFFER (Size > %dMB)", ramLimitMB)
		}
		processOnDisk(ctx, client, bucket, s3Key, startByte, endByte, outputFifoPath, debug)
	}
}

// processInMemory downloads to a byte slice and writes to FIFO
func processInMemory(ctx context.Context, client *s3.Client, bucket, key string, start, end int64, fifoPath string, debug bool) {
	size := end - start + 1
	
	// Allocate buffer
	// We use WriteAtBuffer because S3 Downloader requires io.WriterAt interface for parallel downloads
	buffer := manager.NewWriteAtBuffer(make([]byte, size))

	downloader := manager.NewDownloader(client, func(d *manager.Downloader) {
		d.PartSize = DefaultChunkSize
		d.Concurrency = 10 // Higher concurrency for RAM
	})

	rangeHeader := fmt.Sprintf("bytes=%d-%d", start, end)
	startDown := time.Now()

	// Download to RAM
	_, err := downloader.Download(ctx, buffer, &s3.GetObjectInput{
		Bucket: aws.String(bucket),
		Key:    aws.String(key),
		Range:  aws.String(rangeHeader),
	})
	if err != nil {
		log.Fatalf("RAM Download failed: %v", err)
	}

	if debug {
		dur := time.Since(startDown)
		speed := (float64(size) / 1024 / 1024) / dur.Seconds()
		log.Printf("[GO-RAM] Downloaded %.2f MB in %.2fs (%.2f MB/s)", float64(size)/1024/1024, dur.Seconds(), speed)
	}

	// Write to FIFO
	writeToFifo(fifoPath, buffer.Bytes(), debug)
}

// processOnDisk downloads to a temp file and uses zero-copy splice to FIFO
func processOnDisk(ctx context.Context, client *s3.Client, bucket, key string, start, end int64, fifoPath string, debug bool) {
	tmpFileName := fmt.Sprintf("shard_%d_%d.tmp", os.Getpid(), start)
	tmpFilePath := filepath.Join(os.TempDir(), tmpFileName)

	tmpFile, err := os.Create(tmpFilePath)
	if err != nil {
		log.Fatalf("Failed to create temp file: %v", err)
	}
	defer func() {
		tmpFile.Close()
		os.Remove(tmpFilePath)
	}()

	downloader := manager.NewDownloader(client, func(d *manager.Downloader) {
		d.PartSize = DefaultChunkSize
		d.Concurrency = 5 // Slightly lower concurrency for Disk to avoid IO contention
	})

	rangeHeader := fmt.Sprintf("bytes=%d-%d", start, end)
	startDown := time.Now()

	_, err = downloader.Download(ctx, tmpFile, &s3.GetObjectInput{
		Bucket: aws.String(bucket),
		Key:    aws.String(key),
		Range:  aws.String(rangeHeader),
	})
	if err != nil {
		log.Fatalf("Disk Download failed: %v", err)
	}
	tmpFile.Sync()

	if debug {
		dur := time.Since(startDown)
		log.Printf("[GO-DISK] Download complete in %.2fs", dur.Seconds())
	}

	// Rewind file for reading
	tmpFile.Seek(0, 0)
	
	// Open FIFO and stream
	writeToFifoFromReader(fifoPath, tmpFile, debug)
}

// writeToFifo writes a byte slice to the pipe
func writeToFifo(fifoPath string, data []byte, debug bool) {
	if debug {
		log.Printf("[GO-IO] Opening FIFO for write...")
	}
	f, err := os.OpenFile(fifoPath, os.O_WRONLY, os.ModeNamedPipe)
	if err != nil {
		log.Fatalf("Failed to open FIFO: %v", err)
	}
	defer f.Close()

	if debug {
		log.Printf("[GO-IO] Writing to FIFO...")
	}

	_, err = f.Write(data)
	if err != nil {
		handlePipeError(err, debug)
	}
}

// writeToFifoFromReader uses io.Copy (sendfile)
func writeToFifoFromReader(fifoPath string, src io.Reader, debug bool) {
	if debug {
		log.Printf("[GO-IO] Opening FIFO for write...")
	}
	f, err := os.OpenFile(fifoPath, os.O_WRONLY, os.ModeNamedPipe)
	if err != nil {
		log.Fatalf("Failed to open FIFO: %v", err)
	}
	defer f.Close()

	if debug {
		log.Printf("[GO-IO] Streaming to FIFO (Zero-Copy)...")
	}

	_, err = io.Copy(f, src)
	if err != nil {
		handlePipeError(err, debug)
	}
}

func handlePipeError(err error, debug bool) {
	// Ignore broken pipe (EPIPE) which happens if downstream reader (e.g. 'head') closes early
	if pathErr, ok := err.(*os.PathError); ok && pathErr.Err == syscall.EPIPE {
		if debug {
			log.Printf("[GO-IO] Downstream closed early (EPIPE)")
		}
		return
	}
	if err == syscall.EPIPE {
		return
	}
	log.Fatalf("Write error: %v", err)
}

// Helpers
func parseArgs(args []string) map[string]string {
	m := make(map[string]string)
	for _, arg := range args {
		parts := strings.SplitN(arg, "=", 2)
		if len(parts) == 2 {
			m[parts[0]] = parts[1]
		}
	}
	return m
}

func parseByteRange(s string) (int64, int64, error) {
	if !strings.HasPrefix(s, "bytes=") {
		return 0, 0, fmt.Errorf("must start with 'bytes='")
	}
	parts := strings.Split(s[6:], "-")
	if len(parts) != 2 {
		return 0, 0, fmt.Errorf("invalid format")
	}
	start, err := strconv.ParseInt(parts[0], 10, 64)
	if err != nil {
		return 0, 0, err
	}
	end, err := strconv.ParseInt(parts[1], 10, 64)
	if err != nil {
		return 0, 0, err
	}
	return start, end, nil
}