## Making a new PaSh release

There are two types of PaSh releases. Big feature-full ones that advance PaSh's major version, e.g., `0.7` to `0.8`, and originate from `future`.
There are also small releases that mostly contain bug fixes that advance the minor version,
e.g., `0.7` to `0.7.1`, and go straight to `main`.

For both of them, before doing anything else, we need to modify `compiler/config.py` to update the version. This is necessary because it is parsed to generate the website and other things.

Then we would merge the PR to `main` (either from `future` or from a bug-fix branch).

We then use Github's release feature to create a release, we generate an automatic changelog, and then filter out the important commits/PRs.

This completes the release, and then the CI should:
- Push docker images on docker.io and github images
- Run tests and update badges on main
- Generate and push the new website