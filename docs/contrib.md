## PaSh on Docker: A Pocket Guide

This is a pocket guide for running PaSh in a docker container.

### Loading image

The latest docker container runs PaSh in a pre-configured Ubuntu 18.04 image, loadable as 

```sh
curl img.pash.ndr.md | docker load
```

Then, to create a writable container from this image, run:

```sh
docker run --name pash-playground -it pash/18.04  
```

This will use `pash-playground` for the container (to be able to start/stop it later) and `-it` runs it interactively (subsuming the next command).

After exiting (or if not `-it` is used), one can start the container with:
```sh
docker start -i pash-playground
```

Flag `-i` starts it interactively.

### Customizing image

To be used for continuous integration and testing, this image has been configured to have _read-only_ access to the repo through a different user.
To get write access, update the name and email used for `git` commits with your name and email, and [generate an ssh key pair](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```sh
git config --global user.name "FIRST_NAME LAST_NAME"
git config --global user.email "YOUR_EMAIL@example.com"
ssh-keygen -t rsa -b 4096 -C "YOUR_EMAIL@example.com"
```

Finally, [add the key to your GitHub account](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account).

### Basic Docker commands

This is a small list of useful docker commands. Docker distinguishes between an image (a pre-setup environment) and a container (a running or stopped instance of an image, created with `run`). `NN` below refers to the name of the container, provided to `run`. Containers are always identified uniquely by a sha256 string, which can be useful if no name is given.

```
docker images                    # shows all images in the system
docker ps -a | grep pash         # shows all containers with name pash
docker run --name NN pash/18.04  # creates a wriable container (add `-it` for interactive)
docker start -i NN               # starts container (add `-i` to drop straight into shell)
docker attach NN                 # non-interactive -> interactive
<CTL+p>, then <CTL+q>            # interactive -> non-interactive
docker cp A B                    # copy host<->container; A B can be `NN:/x/y/z` or `lo/c/al`
```

Useful options for `docker run`:
* [Mount host storage](https://docs.docker.com/storage/bind-mounts/): `-v /HST:/IN/NN`
* [Limit CPU/Mem](https://docs.docker.com/config/containers/resource_constraints/): `--cpus='.5' --mem=1g`
* [Limit disk size](https://docs.docker.com/engine/reference/commandline/run/#set-storage-driver-options-per-container): `--storage-opt size=10G`

## Using GNU Screen

We often use GNU `screen`, a terminal multiplexer, for pair programming. GNU `screen` has two benefits: (1) it allows a session to continue executing even after one detaches and closes the terminal window, which is useful for executing long-running processes such as PaSh's experiments; and (2) it allows multiple people to share a session, with the ability to read and write in the terminal window as if it was a single person.

Here are commands related to launching a screen session:
* `screen`                   -> start a new session
* `screen -ls`               -> show all screen sessions in this machine
* `screen -x <ID>`           -> attach to screen with id `<ID>`, as shown by `-ls` above.
* `screen -x <user>/<pash>`  -> attach to session `pash` of user `user`, assuming it exists/running; 

When in a `screen` session, all `screen`-related commands are prefixed by `ctr-a` (which means pressing `ctrl` and `a` _together_, and _then_ pressing the followup character). Here are the 5 most useful commands:
* `ctrl-a c`                 -> create­a new window in the current session
* `ctrl-a "`                 -> show all windows and choose using up/down arrows
* `ctr-a <N>`                -> go to window <N>, where <N> is 0..9 inclusive
* `ctrl-a ctrl-a`            -> fast-toggle between latest windows
* `ctrl-a d`                 -> detach window, to exit without closing screen

A slightly more extended pocket guide [in this gist](https://gist.github.com/nvasilakis/826e4f88d0e0dba2adf4df4834cb9394).

## Pushing to new branches

To keep the main branch free of breakage, we often push technical work in new branches. To do this, from the main/master branch on your local machine:

```sh
git checkout -b <NAME>    # create a new branch and switch to it
git push -u origin <NAME> # push the new branch as a new branch on github
```

You can always fetch changes into the main/master:

```sh
git checkout master        # switch to main/master
git pull                   # fetch remote changes
git checkout <NAME>        # switch back to your branch
git merge master           # fetch changes from main/master
```

(You can use `rebase` instead of `merge` if your branch is local and hasn't been pushed to GitHub, but `merge` if your branch is already pushed.)

