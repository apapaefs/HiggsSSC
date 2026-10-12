# A few terminal commands

Run simulations on Timur inside `screen` so a dropped SSH connection does not
end the shell:

```bash
ssh YOUR_USERNAME@timur.kennesaw.edu
screen -S higgsssc
cd ~/HiggsSSC
```

Run the pipeline there. Detach with **Ctrl-a, then d**; reconnect with
`screen -r higgsssc`. `screen -ls` lists sessions. Do not start the same run in
two sessions.

| Command | Use |
| --- | --- |
| `pwd` | Show the current directory. |
| `ls` | List its files. |
| `cd PATH` | Change directory. |
| `less FILE` | Read a log; press `q` to leave. |
| `tail -n 40 FILE` | Read the most recent log lines. |
| `tail -f FILE` | Follow a log; Ctrl-c stops following, not the simulation. |
| `top -u "$USER"` | Inspect your process activity. |

Use an editor you know (`nano FILE` is a simple option). Preserve the line breaks
in copied command blocks. A trailing backslash continues a shell command;
Python here-documents require their closing marker on its own line.

For reports, see [remote viewing](../SSHFS_REMOTE_MOUNT.md). Avoid deleting
campaign files to recover disk space until you have checked which logs,
manifests and inputs are needed for reproducibility.
