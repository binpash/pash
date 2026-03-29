#!/usr/bin/env bash
set -euo pipefail

IMAGE_BUILDER_DIR="${IMAGE_BUILDER_DIR:-/opt/image-builder}"
BUNDLE_ROOT="${BUNDLE_ROOT:-/opt/minimal-runtime}"
RUNTIME_DIR="${RUNTIME_DIR:-$BUNDLE_ROOT/var/task/runtime}"
PRE_INSTALL_PACKAGES_FILE="${PRE_INSTALL_PACKAGES_FILE:-/tmp/pre-install-packages.txt}"
POST_INSTALL_PACKAGES_FILE="${POST_INSTALL_PACKAGES_FILE:-/tmp/post-install-packages.txt}"
NEW_PACKAGES_FILE="${NEW_PACKAGES_FILE:-/tmp/new-packages.txt}"

mkdir -p "$RUNTIME_DIR"

read_list_file() {
  local file_path=$1
  [[ -f "$file_path" ]] || return 0

  awk '
    {
      sub(/[[:space:]]*#.*/, "", $0)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", $0)
      if (length($0) > 0) {
        print $0
      }
    }
  ' "$file_path"
}

binary_name_from_entry() {
  local entry=$1
  printf '%s\n' "${entry%%=*}"
}

package_name_from_entry() {
  local entry=$1
  if [[ "$entry" == *=* ]]; then
    printf '%s\n' "${entry#*=}"
  else
    printf '%s\n' "$entry"
  fi
}

is_local_runtime_binary() {
  case "$1" in
    split|eager|r_merge|r_split|r_wrap|r_unwrap|set-diff|dgsh-tee)
      return 0
      ;;
    *)
      return 1
      ;;
  esac
}

resolve_package_name() {
  local binary_name=$1
  local package_name

  package_name="$(dnf repoquery --whatprovides "/usr/bin/$binary_name" --qf '%{name}' 2>/dev/null | head -n 1 || true)"
  if [[ -n "$package_name" ]]; then
    printf '%s\n' "$package_name"
    return 0
  fi

  package_name="$(dnf repoquery --whatprovides "*/$binary_name" --qf '%{name}' 2>/dev/null | head -n 1 || true)"
  if [[ -n "$package_name" ]]; then
    printf '%s\n' "$package_name"
    return 0
  fi

  package_name="$(dnf repoquery "$binary_name" 2>/dev/null | sed -n 's/-[0-9].*$//p' | head -n 1 || true)"
  if [[ -n "$package_name" ]]; then
    printf '%s\n' "$package_name"
    return 0
  fi

  return 1
}

resolve_binary_path() {
  local binary_name=$1
  local candidate

  candidate="$(command -v "$binary_name" || true)"
  if [[ -n "$candidate" ]]; then
    printf '%s\n' "$candidate"
    return 0
  fi

  for candidate in \
    "/usr/bin/$binary_name" \
    "/usr/sbin/$binary_name" \
    "/bin/$binary_name" \
    "/sbin/$binary_name"
  do
    if [[ -x "$candidate" ]]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  return 1
}

copy_with_parents() {
  local source_path=$1

  if [[ ! -e "$source_path" && ! -L "$source_path" ]]; then
    return 0
  fi

  (cd / && cp -a --parents ".${source_path}" "$BUNDLE_ROOT")
}

capture_installed_packages() {
  rpm -qa --qf '%{NAME}\n' | sort -u
}

copy_package_payload() {
  local package_name=$1
  local package_path

  while IFS= read -r package_path; do
    [[ -n "$package_path" ]] || continue
    copy_with_parents "$package_path"
  done < <(rpm -ql "$package_name" 2>/dev/null || true)
}

declare -A COPIED_PACKAGE_PAYLOADS=()

copy_package_payload_once() {
  local package_name=$1
  if [[ -n "${COPIED_PACKAGE_PAYLOADS[$package_name]:-}" ]]; then
    return 0
  fi

  copy_package_payload "$package_name"
  COPIED_PACKAGE_PAYLOADS["$package_name"]=1
}

mapfile -t REQUESTED_ENTRIES < <(read_list_file "$IMAGE_BUILDER_DIR/binaries.txt")

SYSTEM_PACKAGES=(binutils)
REQUESTED_PACKAGE_NAMES=()
for entry in "${REQUESTED_ENTRIES[@]}"; do
  binary_name="$(binary_name_from_entry "$entry")"
  if is_local_runtime_binary "$binary_name"; then
    continue
  fi

  requested_package_name="$(package_name_from_entry "$entry")"
  package_name="$(resolve_package_name "$requested_package_name" || true)"
  if [[ -z "$package_name" ]]; then
    echo "Could not resolve an installable package for binary '$binary_name'." >&2
    exit 1
  fi
  SYSTEM_PACKAGES+=("$package_name")
  REQUESTED_PACKAGE_NAMES+=("$package_name")
done

mapfile -t SYSTEM_PACKAGES < <(printf '%s\n' "${SYSTEM_PACKAGES[@]}" | sort -u)
capture_installed_packages > "$PRE_INSTALL_PACKAGES_FILE"
if [[ ${#SYSTEM_PACKAGES[@]} -gt 0 ]]; then
  dnf install -y --setopt=install_weak_deps=0 --nodocs "${SYSTEM_PACKAGES[@]}"
  dnf clean all
fi
capture_installed_packages > "$POST_INSTALL_PACKAGES_FILE"
comm -13 "$PRE_INSTALL_PACKAGES_FILE" "$POST_INSTALL_PACKAGES_FILE" > "$NEW_PACKAGES_FILE"

while IFS= read -r package_name; do
  [[ -n "$package_name" ]] || continue
  copy_package_payload_once "$package_name"
done < "$NEW_PACKAGES_FILE"

for package_name in "${REQUESTED_PACKAGE_NAMES[@]}"; do
  copy_package_payload_once "$package_name"
done

for entry in "${REQUESTED_ENTRIES[@]}"; do
  binary_name="$(binary_name_from_entry "$entry")"
  binary_path="$(resolve_binary_path "$binary_name" || true)"
  if [[ -z "$binary_path" ]]; then
    echo "Binary '$binary_name' is not available after installation." >&2
    exit 1
  fi

  ln -sf "$binary_path" "$RUNTIME_DIR/$binary_name"
done
