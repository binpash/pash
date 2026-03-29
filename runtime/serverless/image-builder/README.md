# Lambda Image Builder

This directory builds an AWS Lambda container image instead of a zip bundle.

The image approach fixes the main package-churn problem from the zip workflow:
- each build starts from a clean base image
- removed Python or system packages do not linger in the artifact
- local runtime binaries are rebuilt from source during the image build

## Files

- `Dockerfile`: multi-stage Lambda image build
- `build-image.sh`: stages a minimal Docker build context and builds the image
- `binaries.txt`: package-manager binaries, one per line
- `python-packages.txt`: pip packages, one per line
- `stage-runtime-assets.sh`: builder-stage script that compiles local runtime binaries from source
- `install-system-binaries.sh`: packager-stage script that installs non-local binaries in a temporary image, copies the full payload of the newly installed RPMs into a bundle, and symlinks the requested commands into `/var/task/runtime`

## What gets built

The image always copies the top-level runtime source files from [`/home/ubuntu/pash/runtime`](/home/ubuntu/pash/runtime) into `/var/task/runtime`:
- `*.c`
- `*.h`
- `*.sh`
- `Makefile`

After copying those files, the image build runs `make` in `/var/task/runtime`, so local runtime binaries such as `r_split`, `r_wrap`, `r_unwrap`, `r_merge`, `split`, `eager`, and `dgsh-tee` are built from source automatically.

For entries in [`binaries.txt`](/home/ubuntu/pash/runtime/serverless/image-builder/binaries.txt):
- non-local binaries such as `openssl` are resolved from Amazon Linux packages in a temporary packager stage
- the full installed payload of those newly added RPMs is copied into the final image, so packages such as `ImageMagick` keep all of their runtime files together
- each requested binary is then symlinked into `/var/task/runtime` so existing scripts can keep using the runtime path

For entries in `python-packages.txt`:
- packages are installed into the Lambda Python environment with `pip`

The Lambda app code copied into the image is:
- [`lambda-function.py`](/home/ubuntu/pash/runtime/serverless/lambda-function.py)
- [`test.sh`](/home/ubuntu/pash/runtime/serverless/test.sh)
- [`aws/`](/home/ubuntu/pash/runtime/serverless/aws)

Those Lambda app files are copied at the end of the final image stage so edits to Python handler code or `test.sh` do not invalidate the expensive runtime build, RPM payload bundling, or pip-install layers.

## Build

```bash
sudo /home/ubuntu/pash/runtime/serverless/image-builder/build-image.sh pash-serverless:latest
```

## Deploy

`setup-lambda-image.sh` uses your existing shell environment directly. Export the required values once, then run the script without an inline env prefix:

```bash
export AWS_ACCOUNT_ID=...
export AWS_BUCKET=...
export AWS_REGION=us-east-1

/home/ubuntu/pash/runtime/serverless/image-builder/setup-lambda-image.sh
```

`REGION` is also accepted, but if it is unset the script now falls back to `AWS_REGION`, then `AWS_DEFAULT_REGION`, then `us-east-1`.

## Runtime layout

The final image keeps runtime helpers and runtime source files under `/var/task/runtime` so existing scripts can keep using the same relative paths.

The image sets:

```bash
PATH=/var/task/runtime:$PATH
```

Unlike the zip workflow, there is no separate `runtime/lib` or `runtime/python` payload to maintain.

## Notes

- `dgsh-tee` is still built through [`/home/ubuntu/pash/runtime/Makefile`](/home/ubuntu/pash/runtime/Makefile), which clones the upstream `dgsh` source during the build.
- If you remove a package from `python-packages.txt` or `binaries.txt` and rebuild the image, it disappears from the new artifact automatically.
- `convert=ImageMagick` maps the runtime command name (`convert`) to the Amazon Linux package name (`ImageMagick`). `convert` is the binary you want on `PATH`; `ImageMagick` is the RPM that provides it.
- The final image should not contain `dnf`, `gcc`, `git`, or `make`; those remain confined to build-time stages.
