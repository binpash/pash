# Stage 1: build the local runtime binaries in a toolchain image. This stage
# is intentionally heavy because it needs gcc/git/make, but none of that should
# ship in the final Lambda image.
FROM public.ecr.aws/amazonlinux/amazonlinux:2023 AS runtime-builder

RUN dnf install -y \
      gcc \
      gcc-c++ \
      git \
      make \
      cargo \
      findutils \
    && dnf clean all

WORKDIR /src

COPY runtime /src/runtime
COPY serverless-runtime /src/serverless-runtime
COPY serverless-pashlib /src/serverless-pashlib
COPY image-builder /src/image-builder

RUN /bin/bash /src/image-builder/stage-runtime-assets.sh

# Stage 2: install package-manager binaries in a temporary image and copy the
# full payload of the newly installed RPMs into a bundle. This keeps dnf and
# package-manager state out of the final image while preserving complete
# runtime payloads for packages such as ImageMagick.
FROM public.ecr.aws/lambda/python:3.12 AS system-binary-packager

RUN dnf install -y \
      binutils \
      findutils \
    && dnf clean all

WORKDIR /opt

COPY image-builder /opt/image-builder

RUN /bin/bash /opt/image-builder/install-system-binaries.sh

# Stage 3: start again from a clean Lambda Python base image and copy in only
# the artifacts we want to deploy. Multi-stage builds do not share files
# automatically, so COPY --from is what keeps the final image small and avoids
# shipping compilers, package-manager state, and other build-time baggage.
FROM public.ecr.aws/lambda/python:3.12

WORKDIR /var/task

COPY image-builder/python-packages.txt /opt/image-builder/python-packages.txt
COPY image-builder/binaries.txt /opt/image-builder/binaries.txt
COPY image-builder/verify-runtime-linking.sh /opt/image-builder/verify-runtime-linking.sh
COPY vendored-runtime /opt/vendored-runtime
COPY --from=runtime-builder /opt/pash-built/runtime /var/task/runtime
COPY --from=system-binary-packager /opt/minimal-runtime/ /

RUN cp -Lf /opt/vendored-runtime/* /var/task/runtime/ \
    && if grep -Eq '^[[:space:]]*[^#[:space:]]' /opt/image-builder/python-packages.txt 2>/dev/null; then \
         python3 -m pip install --no-cache-dir -r /opt/image-builder/python-packages.txt; \
       fi \
    && /bin/bash -lc ' \
         shopt -s nullglob; \
         for file in /var/task/runtime/*.sh; do \
           chmod +x "$file"; \
         done; \
         for file in /var/task/runtime/*; do \
           [[ -f "$file" ]] || continue; \
           case "$file" in \
             *.py|*.sh) ;; \
             *) chmod +x "$file" ;; \
           esac; \
         done' \
    && rm -rf /root/.cache/pip

COPY serverless/aws /var/task/aws
COPY serverless/lambda-function.py /var/task/lambda-function.py
COPY serverless/test.sh /var/task/test.sh

RUN chmod +x /var/task/test.sh

ENV PATH="/var/task/runtime:${PATH}"
ENV LD_LIBRARY_PATH="/var/task/runtime/lib:${LD_LIBRARY_PATH}"

CMD ["lambda-function.lambda_handler"]
