package main

import (
	"context"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"
	"syscall"
	"time"

	"github.com/aws/aws-sdk-go-v2/aws"
	"github.com/aws/aws-sdk-go-v2/config"
	"github.com/aws/aws-sdk-go-v2/service/s3"
)

// Global constants
const (
	// Use a 1MB buffer for the io.Copy operation to maximize throughput
	// Default io.Copy uses 32KB, which might be too small for high-bandwidth S3 streaming
	CopyBufferSize = 1 * 1024 * 1024
)

// Timing infrastructure for performance debugging
type timingEvent struct {
	stage       string
	elapsedMs   int
	description string
}

var timingStart = time.Now()
var timingLog = make([]timingEvent, 0)

func logTiming(stage, description string, debug bool) {
	elapsed := time.Since(timingStart)
	event := timingEvent{
		stage:       stage,
		elapsedMs:   int(elapsed.Milliseconds()),
		description: description,
	}
	timingLog = append(timingLog, event)
	if debug {
		log.Printf("[%.3fs][TIMING] %s: %s", elapsed.Seconds(), stage, description)
	}
}

func printTimingSummary(debug bool) {
	if debug && len(timingLog) > 0 {
		log.Println("======================================================================")
		log.Println("[TIMING SUMMARY]")
		log.Println("======================================================================")
		prevMs := 0
		for _, event := range timingLog {
			delta := event.elapsedMs - prevMs
			log.Printf("  %-20s @ %6dms (+%5dms) %s", event.stage, event.elapsedMs, delta, event.description)
			prevMs = event.elapsedMs
		}
		log.Println("======================================================================")
	}
}

func main() {
	// 1. Argument Parsing matching Python signature
	// Usage: <s3_key> <output_fifo> <byte_range> [split=N] [num_shards=N] [job_uid=UID]
	if len(os.Args) < 4 {
		log.Fatalf("Usage: %s <s3_key> <output_fifo> <byte_range> [key=val...]", os.Args[0])
	}

	s3Key := os.Args[1]
	outputFifoPath := os.Args[2]
	byteRangeStr := os.Args[3]

	// Parse keyword args
	kvArgs := make(map[string]string)
	for _, arg := range os.Args[4:] {
		parts := strings.SplitN(arg, "=", 2)
		if len(parts) == 2 {
			kvArgs[parts[0]] = parts[1]
		}
	}

	bucket := os.Getenv("AWS_BUCKET")
	if bucket == "" {
		log.Fatal("Error: AWS_BUCKET environment variable not set")
	}

	// PERF_MODE disables debug output for clean performance tests
	perfMode := os.Getenv("PASH_PERF_MODE") == "true"
	debug := kvArgs["debug"] == "true" && !perfMode

	logTiming("INIT", fmt.Sprintf("Go shard reader started: %s", s3Key), debug)

	if debug {
		log.Printf("[GO-MAIN] S3: %s/%s Range: %s", bucket, s3Key, byteRangeStr)
	}

	// 2. Parse Byte Range
	startByte, endByte, err := parseByteRange(byteRangeStr)
	if err != nil {
		log.Fatalf("Invalid byte range: %v", err)
	}
	// contentLength := endByte - startByte + 1 // Unused in streaming mode

	// 3. Initialize AWS Session
	ctx := context.TODO()
	cfg, err := config.LoadDefaultConfig(ctx)
	if err != nil {
		log.Fatalf("unable to load SDK config, %v", err)
	}
	client := s3.NewFromConfig(cfg)

	// 4. Open FIFO
	// We open this FIRST. This blocks until the downstream reader connects.
	// This mirrors the Python "Direct Pipe" strategy.
	logTiming("FIFO_OPEN_START", fmt.Sprintf("Opening FIFO %s", outputFifoPath), debug)

	if debug {
		log.Printf("[GO-STREAM] Opening FIFO for write (blocking): %s", outputFifoPath)
	}

	fifoFile, err := os.OpenFile(outputFifoPath, os.O_WRONLY, os.ModeNamedPipe)
	if err != nil {
		log.Fatalf("Failed to open FIFO: %v", err)
	}
	defer fifoFile.Close()

	logTiming("FIFO_OPEN_END", "FIFO connected", debug)

	// 5. Direct Streaming: S3 -> FIFO
	// We use GetObject with Range, then CopyBuffer to stream the response body directly.
	logTiming("STREAM_START", "Starting S3 -> FIFO stream", debug)

	if debug {
		log.Printf("[GO-STREAM] Strategy: DIRECT PIPE (Low Latency)")
	}

	rangeHeader := fmt.Sprintf("bytes=%d-%d", startByte, endByte)

	// Execute the request (Headers only received initially)
	resp, err := client.GetObject(ctx, &s3.GetObjectInput{
		Bucket: aws.String(bucket),
		Key:    aws.String(s3Key),
		Range:  aws.String(rangeHeader),
	})
	if err != nil {
		log.Fatalf("Failed to start S3 download: %v", err)
	}
	defer resp.Body.Close()

	startStream := time.Now()

	// 6. IO Copy with explicit buffer
	// This reads chunks from the network socket (resp.Body) and writes them to the FIFO (fifoFile)
	// No intermediate file is used. Memory usage is limited to the buffer size.
	buf := make([]byte, CopyBufferSize)
	written, err := io.CopyBuffer(fifoFile, resp.Body, buf)

	if err != nil {
		// Broken pipe is common if downstream reader exits early (e.g., 'head')
		if pathErr, ok := err.(*os.PathError); ok && pathErr.Err == syscall.EPIPE {
			if debug {
				log.Printf("[GO-STREAM] Downstream closed pipe early (EPIPE)")
			}
		} else {
			log.Fatalf("Failed to write to FIFO: %v", err)
		}
	}

	streamDur := time.Since(startStream)
	speedMB := (float64(written) / 1024 / 1024) / streamDur.Seconds()

	logTiming("STREAM_END", fmt.Sprintf("%d bytes in %.2fs (%.2f MB/s)", written, streamDur.Seconds(), speedMB), debug)

	if debug {
		log.Printf("[GO-STREAM] Streamed %d bytes to FIFO in %.2fs", written, streamDur.Seconds())
	}

	logTiming("COMPLETE", "Go shard reader complete", debug)
	printTimingSummary(debug)
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
