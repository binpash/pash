package main

import (
	"context"
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"runtime/debug"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"
	"github.com/aws/aws-sdk-go-v2/feature/s3/manager"
	"github.com/aws/aws-sdk-go-v2/service/s3"
)

// Global constants
const (
	// Download Chunk Size (16MB is a sweet spot for high-bandwidth links like Lambda/EC2)
	DefaultPartSize = 16 * 1024 * 1024

	// F_SETPIPE_SZ is the fcntl command to change pipe buffer size on Linux
	F_SETPIPE_SZ = 1031
)

// --- Timing Helpers ---
type timingEvent struct {
	stage       string
	elapsedMs   int
	description string
}

var timingStart = time.Now()
var timingLog = make([]timingEvent, 0)

func logTiming(stage, description string, enabled bool) {
	if !enabled {
		return
	}
	elapsed := time.Since(timingStart)
	event := timingEvent{
		stage:       stage,
		elapsedMs:   int(elapsed.Milliseconds()),
		description: description,
	}
	timingLog = append(timingLog, event)
}

func printTimingSummary(enabled bool) {
	if !enabled || len(timingLog) == 0 {
		return
	}
	log.Println("[TIMING SUMMARY]")
	prevMs := 0
	for _, event := range timingLog {
		delta := event.elapsedMs - prevMs
		log.Printf("  %-20s @ %6dms (+%5dms) %s", event.stage, event.elapsedMs, delta, event.description)
		prevMs = event.elapsedMs
	}
}

func main() {
	// 0. OPTIMIZATION: Disable GC initially
	// This prevents pauses during startup and import.
	oldGC := debug.SetGCPercent(-1)

	// 1. Argument Parsing
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

	// Performance/Debug flags
	perfMode := os.Getenv("PASH_PERF_MODE") == "true"
	isDebug := kvArgs["debug"] == "true" && !perfMode

	logTiming("INIT", fmt.Sprintf("Started: %s", s3Key), isDebug)

	// Determine RAM buffer threshold (default 128MB)
	ramLimitMB := 128
	if envLimit := os.Getenv("PASH_S3_RAM_BUFFER_MB"); envLimit != "" {
		if val, err := strconv.Atoi(envLimit); err == nil {
			ramLimitMB = val
		}
	}
	ramLimitBytes := int64(ramLimitMB) * 1024 * 1024

	// 2. Parse Byte Range
	startByte, endByte, err := parseByteRange(byteRangeStr)
	if err != nil {
		log.Fatalf("Invalid byte range: %v", err)
	}
	contentLength := endByte - startByte + 1

	// 3. Initialize AWS Session with Transport Tuning
	// Custom HTTP client for raw throughput
	httpClient := &http.Client{
		Transport: &http.Transport{
			DisableCompression:  true, // Save CPU
			MaxIdleConns:        100,
			MaxIdleConnsPerHost: 100,
			IdleConnTimeout:     90 * time.Second,
		},
	}

	ctx := context.TODO()
	cfg, err := config.LoadDefaultConfig(ctx, config.WithHTTPClient(httpClient))
	if err != nil {
		log.Fatalf("unable to load SDK config, %v", err)
	}
	client := s3.NewFromConfig(cfg)
	logTiming("AWS_INIT", "SDK Loaded", isDebug)

	// 4. Execution Strategy
	if contentLength <= ramLimitBytes {
		// Strategy A: IN-MEMORY (Fastest)
		if isDebug {
			log.Printf("[GO-FAST] Mode: RAM (Size %.2f MB <= %d MB limit)", float64(contentLength)/1024/1024, ramLimitMB)
		}
		// Pass oldGC to restore if needed (though we keep it off for speed)
		processInMemory(ctx, client, bucket, s3Key, startByte, endByte, outputFifoPath, isDebug, oldGC)
	} else {
		// Strategy B: DISK BUFFER (Safe)
		if isDebug {
			log.Printf("[GO-FAST] Mode: DISK (Size %.2f MB > %d MB limit)", float64(contentLength)/1024/1024, ramLimitMB)
		}
		// Re-enable GC for disk mode to prevent OOM on long running streams,
		// but keep it conservative (100% is default)
		debug.SetGCPercent(100)
		processOnDisk(ctx, client, bucket, s3Key, startByte, endByte, outputFifoPath, isDebug)
	}

	printTimingSummary(isDebug)
}

// processInMemory downloads to RAM buffer -> Writes to FIFO
func processInMemory(ctx context.Context, client *s3.Client, bucket, key string, start, end int64, fifoPath string, isDebug bool, oldGC int) {
	size := end - start + 1

	// Pre-allocate buffer. Since GC is off, this is just a big malloc.
	data := make([]byte, size)
	buffer := manager.NewWriteAtBuffer(data)

	// Setup parallel downloader
	downloader := manager.NewDownloader(client, func(d *manager.Downloader) {
		d.PartSize = DefaultPartSize
		d.Concurrency = 15 // High concurrency for RAM
		// Pooled buffer provider reduces allocation churn
		d.BufferProvider = manager.NewPooledBufferedWriterReadFromProvider(int(DefaultPartSize))
	})

	rangeHeader := fmt.Sprintf("bytes=%d-%d", start, end)

	logTiming("DL_START", "Starting RAM download", isDebug)

	_, err := downloader.Download(ctx, buffer, &s3.GetObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key), Range: aws.String(rangeHeader),
	})
	if err != nil {
		log.Fatalf("RAM Download failed: %v", err)
	}

	logTiming("DL_END", fmt.Sprintf("Downloaded %d bytes to RAM", size), isDebug)

	// Open FIFO and Write
	writeToFifoWithTuning(fifoPath, buffer.Bytes(), isDebug)
}

// processOnDisk downloads to /tmp -> Zero-copy to FIFO
func processOnDisk(ctx context.Context, client *s3.Client, bucket, key string, start, end int64, fifoPath string, isDebug bool) {
	// Create temp file
	tmpFileName := fmt.Sprintf("shard_%d_%d.tmp", os.Getpid(), start)
	tmpFilePath := filepath.Join(os.TempDir(), tmpFileName)

	tmpFile, err := os.Create(tmpFilePath)
	if err != nil {
		log.Fatalf("Tmp create failed: %v", err)
	}

	// Cleanup ensures disk is cleared
	defer func() {
		tmpFile.Close()
		os.Remove(tmpFilePath)
	}()

	downloader := manager.NewDownloader(client, func(d *manager.Downloader) {
		d.PartSize = DefaultPartSize
		d.Concurrency = 5 // Lower concurrency for disk to avoid IO contention
	})

	logTiming("DL_START", "Starting Disk download", isDebug)

	_, err = downloader.Download(ctx, tmpFile, &s3.GetObjectInput{
		Bucket: aws.String(bucket), Key: aws.String(key), Range: aws.String(fmt.Sprintf("bytes=%d-%d", start, end)),
	})
	if err != nil {
		log.Fatalf("Disk Download failed: %v", err)
	}

	// Ensure data is on disk
	tmpFile.Sync()
	logTiming("DL_END", "Download complete", isDebug)

	// Rewind for reading
	tmpFile.Seek(0, 0)

	// Stream using zero-copy
	writeToFifoFromReaderWithTuning(fifoPath, tmpFile, isDebug)
}

// writeToFifoWithTuning handles opening FIFO, tuning kernel buffer, and writing data
func writeToFifoWithTuning(fifoPath string, data []byte, isDebug bool) {
	logTiming("FIFO_OPEN", "Opening FIFO...", isDebug)

	f, err := os.OpenFile(fifoPath, os.O_WRONLY, os.ModeNamedPipe)
	if err != nil {
		log.Fatalf("FIFO open failed: %v", err)
	}
	defer f.Close()

	// Tune Kernel Pipe Buffer
	tunePipeBuffer(f, isDebug)

	logTiming("WRITE_START", "Flushing RAM to FIFO", isDebug)
	_, err = f.Write(data)
	if err != nil {
		handlePipeError(err, isDebug)
	}
	logTiming("WRITE_END", "Flush complete", isDebug)
}

// writeToFifoFromReaderWithTuning handles zero-copy streaming
func writeToFifoFromReaderWithTuning(fifoPath string, src io.Reader, isDebug bool) {
	logTiming("FIFO_OPEN", "Opening FIFO...", isDebug)

	f, err := os.OpenFile(fifoPath, os.O_WRONLY, os.ModeNamedPipe)
	if err != nil {
		log.Fatalf("FIFO open failed: %v", err)
	}
	defer f.Close()

	tunePipeBuffer(f, isDebug)

	logTiming("STREAM_START", "Streaming Disk -> FIFO", isDebug)
	// io.Copy uses splice/sendfile on Linux automatically
	_, err = io.Copy(f, src)
	if err != nil {
		handlePipeError(err, isDebug)
	}
	logTiming("STREAM_END", "Stream complete", isDebug)
}

// tunePipeBuffer attempts to increase the Linux pipe buffer size to 1MB
func tunePipeBuffer(f *os.File, isDebug bool) {
	fd := int(f.Fd())
	// 1048576 = 1MB
	_, _, errno := syscall.Syscall(syscall.SYS_FCNTL, uintptr(fd), F_SETPIPE_SZ, uintptr(1048576))
	if errno != 0 {
		if isDebug {
			log.Printf("[GO-TUNE] Warn: Failed to set pipe buffer: %d", errno)
		}
	} else {
		if isDebug {
			log.Printf("[GO-TUNE] Pipe buffer set to 1MB")
		}
	}
}

func handlePipeError(err error, isDebug bool) {
	// Ignore EPIPE (Broken Pipe) - happens if downstream tool (e.g. head) exits early
	if pathErr, ok := err.(*os.PathError); ok && pathErr.Err == syscall.EPIPE {
		return
	}
	if err == syscall.EPIPE {
		return
	}
	log.Fatalf("Write error: %v", err)
}

// --- Helpers ---

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
		return 0, 0, fmt.Errorf("bad range")
	}
	parts := strings.Split(s[6:], "-")
	start, _ := strconv.ParseInt(parts[0], 10, 64)
	end, _ := strconv.ParseInt(parts[1], 10, 64)
	return start, end, nil
}
