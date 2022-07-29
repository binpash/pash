# Packaging PaSh

`scripts/package` packages PaSh into various formats. The software inside is a suitable
replacement for deployments depending on `install.sh`,
`setup-pash.sh`, `up.sh`, `distro-deps.sh`, and/or `pkg.sh`.


## Building Packages
[eef]: https://github.com/docker/cli/blob/master/experimental/README.md

```
build.sh VERSION [OUTPUT_FORMAT ...]
```

Creates new package files using the repository's current content, each
with mode 440. All output files appear in `scripts/package/output`.

**Clean up junk files and new test results!** The scripts use a
reasonable ignore list, but otherwise defaults to including files
within the PaSh repository.

Command-line arguments correspond to those found in [FPM's
CLI](https://fpm.readthedocs.io/en/v1.14.2/cli-reference.html):

  - `VERSION` is the value of `--version`.
  - `OUTPUT_FORMAT` is some value for `--output-type`


## Testing Packages

```
  deploy.sh IMAGE VERSION FORMAT
```

`deploy.sh` runs `build.sh VERSION FORMAT`, then installs PaSh in a
Docker container based on a `IMAGE` from DockerHub.

In this example, `deploy.sh` installs a `pacman`/`.tar.gz` package
deploys in Arch, a `.deb` package for Debian and Ubuntu, and an RPM
package for Fedora.

```
$ version=0.0.1
$ PATH="$PASH_TOP/scripts/package:$PATH"
$ deploy.sh archlinux $version pacman
$ deploy.sh debian $version deb
$ deploy.sh ubuntu $version deb
$ deploy.sh fedora $version rpm
```

## Developing Packages

```
  iterate.sh IMAGE FORMAT
```

Runs `deploy.sh IMAGE 0.0.1 FORMAT`, then prompts to repeat.


## GitHub Actions

The source code repository defines workflows for GitHub actions.  The
workflows for packaging all use `deploy.sh` to build and run PaSh's
"Hello, World" example. One workflow responds to new pull requests.

**On-Demand Packaging** is the only workflow that can be directly
used in the GitHub GUI. You can find it under the **Actions** tab.

![image](https://user-images.githubusercontent.com/1312121/181671950-89ec5f57-5b9f-4fdb-90a2-6d1099257ad1.png)

Click "Run Workflow" and fill out the form. The form defines what commit of PaSh gets packaged, and each package format may have its own version string.
