# PaSh Troubleshooting

This is a list of issues when setting up or running PaSh (and how to solve them).

### Module `pexpect` not found OR libdash library for parser not working  

_Solution_: run dependency scipt with **sudo** and build pash libraries 

```console
cd $PASH_TOP  
sudo bash scripts/distro-deps.sh  
bash scripts/setup-pash.sh 
```

