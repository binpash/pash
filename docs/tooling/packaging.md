# Packaging PaSh

`scripts/package` packages PaSh into various formats. It is a suitable
replacement for deployments depending on `install.sh`,
`setup-pash.sh`, `up.sh`, `distro-deps.sh`, and/or `pkg.sh`.



## Building Packages
[eef]: https://github.com/docker/cli/blob/master/experimental/README.md

```
  shell.sh
```

To use this feature, [enable experimental features][eef] in Docker.

`shell.sh` enters a Docker container for building packages. Run
without arguments to enter a REPL and print usage information for the
`build` command.

Since `shell.sh` is just a Bash, you can pass it commands.  For
example, `./shell.sh -ic 'build 0.0.1 deb'` launches an interactive
shell that builds a Debian package using version `0.0.1`.


## Testing Packages

```
  deploy.sh IMAGE VERSION FORMAT
```

`deploy.sh` runs `./shell.sh -ic 'build VERSION FORMAT'`, then
installs PaSh in a Docker container based on a `IMAGE` from DockerHub.

Here's an example that builds a `pacman`/`.tar.gz` package for
deployment in Arch, a `.deb` package for Debian and Ubuntu, and an RPM
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

Run `deploy.sh IMAGE 0.0.1 FORMAT`, then prompt to repeat.
