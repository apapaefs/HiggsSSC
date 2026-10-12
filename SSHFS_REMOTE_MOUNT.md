# View reports from Timur

Reports are static HTML with neighboring images and tables. Copy the **whole
report directory**, then open its `index.html` locally:

```bash
scp -r YOUR_USERNAME@timur.kennesaw.edu:~/HiggsSSC/hgammagamma/LOAnalysis/plots/lo_pilot_01 ./lo_pilot_01
```

Use your actual remote checkout path and run tag. For four leptons, copy
`hfourlepton/runs/TAG/report/`. Keep original campaign inputs on Timur;
copying a report does not create a portable analysis campaign.

## Optional Linux SSHFS mount

Install SSHFS using your distribution's package manager (`sshfs` package),
then mount your own checkout read-only:

```bash
mkdir -p ~/mnt/timur-higgsssc
sshfs YOUR_USERNAME@timur.kennesaw.edu:~/HiggsSSC ~/mnt/timur-higgsssc \
  -o ro,reconnect,ServerAliveInterval=15,ServerAliveCountMax=3
xdg-open ~/mnt/timur-higgsssc/hgammagamma/LOAnalysis/plots/lo_pilot_01/index.html
```

Run generators and analysis directly on Timur. The mount is for viewing files.
When finished (or if the mount becomes stale):

```bash
fusermount3 -u ~/mnt/timur-higgsssc
```

Older distributions may use `fusermount -u` instead.
