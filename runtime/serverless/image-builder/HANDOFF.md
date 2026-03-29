# Handoff

## Goal

Move the Lambda packaging flow from zip-based `runtime/` bundling to a container-image workflow.

The desired model is:
- copy top-level runtime source files from [`/home/ubuntu/pash/runtime`](/home/ubuntu/pash/runtime) into the image
- run `make` there to build local runtime binaries
- use [`binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt) only for binaries that should come from the package manager
- use [`python-packages.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/python-packages.txt) for Python dependencies

## Current Files

Image-builder files are under:
- [`/home/ubuntu/pash/runtime/serverless/image-builder/Dockerfile`](/home/ubuntu/pash/runtime/serverless/image-builder/Dockerfile)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/build-image.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/build-image.sh)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/install-system-binaries.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/install-system-binaries.sh)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/stage-runtime-assets.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/stage-runtime-assets.sh)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/python-packages.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/python-packages.txt)
- [`/home/ubuntu/pash/runtime/serverless/image-builder/README.md`](/home/ubuntu/pash/runtime/serverless/image-builder/README.md)

## What Was Changed

- The image builder now copies top-level `*.c`, `*.h`, `*.sh`, and `Makefile` from [`/home/ubuntu/pash/runtime`](/home/ubuntu/pash/runtime) into `/var/task/runtime`.
- [`stage-runtime-assets.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/stage-runtime-assets.sh) now just runs `make` after copying those files.
- [`binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt) is reserved for package-manager binaries only.
- [`README.md`](/home/ubuntu/pash/runtime/serverless/image-builder/README.md) was updated to match that design.
- [`/home/ubuntu/pash/runtime/serverless/README.md`](/home/ubuntu/pash/runtime/serverless/README.md) now points to the image-builder docs.
- [`build-image.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/build-image.sh) now stages its Docker build context under `runtime/serverless/` instead of `/tmp` to avoid the local Docker confinement issue on this machine.
- [`install-system-binaries.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/install-system-binaries.sh) was simplified to install package-manager binaries directly into the final Lambda image and now supports `binary=package` entries in [`binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt).
- [`Dockerfile`](/home/ubuntu/pash/runtime/serverless/image-builder/Dockerfile) now installs package-manager binaries directly in the final image instead of copying a curated minimal runtime subset, and still copies vendored runtime binaries from `/opt/vendored-runtime` into `/var/task/runtime`.
- [`setup-lambda-image.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/setup-lambda-image.sh) was added as a separate image-based deployment script. It does not replace the existing zip-based [`/home/ubuntu/pash/runtime/serverless/setup-lambda.sh`](/home/ubuntu/pash/runtime/serverless/setup-lambda.sh).
- [`test.sh`](/home/ubuntu/pash/runtime/serverless/test.sh) now verifies the installed runtime toolchain instead of only printing the OpenSSL version.
- [`lambda-function.py`](/home/ubuntu/pash/runtime/serverless/lambda-function.py) now runs `test.sh` with `check=True` so test failures propagate to Lambda.

## Current Lists

[`binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt):

```text
openssl
gzip
convert=ImageMagick
file
xargs=findutils
tcpdump
col=util-linux
rev=util-linux
which
```

[`python-packages.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/python-packages.txt):

```text
tensorflow
```

Vendored runtime binaries copied into the image outside the package manager path:

```text
ffmpeg
```

## Validation Done

- Shell syntax checks passed for:
  - [`build-image.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/build-image.sh)
  - [`stage-runtime-assets.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/stage-runtime-assets.sh)
  - [`install-system-binaries.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/install-system-binaries.sh)
- Shell syntax checks passed for:
  - [`setup-lambda-image.sh`](/home/ubuntu/pash/runtime/serverless/image-builder/setup-lambda-image.sh)
- The full Docker image build completed successfully.
- The rebuilt local image size after adding the extra runtime tools is `3,015,943,198` bytes (`~3.02 GB`, `~2.81 GiB`).
- A new image-based Lambda function was created:
  - function name: `lambda-image`
  - region: `us-east-1`
  - package type: `Image`
- `lambda-image` was updated to use the latest ECR image digest:
  - `sha256:640457d6a6b286cb1d2c1a77a70dd637485590a39bc0a715aa1a3928879cd224`
- `lambda-image` was invoked successfully with timeout set to `10` seconds, and the runtime test completed inside Lambda.
- The Lambda-side test verified:
  - `openssl`
  - `ffmpeg`
  - `gzip`
  - `convert`
  - `file`
  - `xargs`
  - `tcpdump`
  - `col`
  - `rev`
  - `which`
  - `eager`
  - `split`
  - `r_merge`
  - `r_split`
  - `r_wrap`
  - `r_unwrap`
  - `set-diff`
  - `dgsh-tee`

## Not Yet Done

- The image-based deployment script currently pushes the image and then immediately attempts both `update-function-code` and `update-function-configuration`. On existing image functions, AWS may still report `ResourceConflictException` on the configuration step because the image update is still in progress. The safe workaround is to wait for `aws lambda wait function-updated` and then apply the configuration update separately.
- The current runtime test is close to the 10-second timeout in Lambda. Keep a larger timeout for normal operation.
- The image is now large enough that cold-start pull time remains a real concern, especially because `convert=ImageMagick` pulls in a large dependency closure.
- The image no longer uses a curated minimal RPM runtime subset, because tools like ImageMagick were too fragile under partial copying.
- The main zip-based Lambda function named `lambda` still exists and remains `PackageType=Zip`.

## Next Step

Run:

```bash
/home/ubuntu/pash/runtime/serverless/image-builder/setup-lambda-image.sh
```

The script uses the existing exported shell environment directly. `REGION` may be set explicitly, but if it is unset the script falls back to `AWS_REGION`, then `AWS_DEFAULT_REGION`, then `us-east-1`.

If you need to verify the deployed image function after a rebuild, use `lambda-image` and wait for the function update to settle before invoking it.

## Expected Risk Areas

- [`/home/ubuntu/pash/runtime/Makefile`](/home/ubuntu/pash/runtime/Makefile) builds `dgsh-tee` by cloning an external repo during `make`
- `ffmpeg` is not available from the default Lambda-image package repositories and is currently copied from the existing vendored binary at [`/home/ubuntu/pash/runtime/serverless/runtime/ffmpeg`](/home/ubuntu/pash/runtime/serverless/runtime/ffmpeg)
- `convert=ImageMagick` materially increases image size and package count
- the image-based deployment script may need a cleaner two-phase update flow for existing image functions to avoid transient AWS update conflicts
