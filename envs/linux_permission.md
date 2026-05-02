# Linux permissions

## Change permissions
- `755` for directories, `644` for files: everyone can read/copy; only owner can modify.
- `750` for directories, `640` for files: only owner and group can access; others cannot.

```bash
# Share with everyone
find folder -type d -exec chmod 755 {} \;
find folder -type f -exec chmod 644 {} \;

# Share with group only
find folder -type d -exec chmod 750 {} \;
find folder -type f -exec chmod 640 {} \;
```

- `x` on a directory = can enter it.
- `r` on a file = can open/copy it.
- Avoid `chmod -R 755 folder` because it makes all files executable.


## Check permissions:

```bash
ls -ld folder            # show symbolic permission, e.g. drwxr-xr-x
stat -c '%a %n' folder   # show numeric permission, e.g. 755
```