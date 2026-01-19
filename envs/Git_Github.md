## Git - Github, Connection

**The Problem**

GitHub authentication fails in every new terminal window unless you manually re-run the setup commands (Using Option 1).

**The Reason**

SSH only looks for standard filenames (like id_rsa) by default. Since your key uses a custom name (id_tiger), the system ignores it and "forgets" it the moment you close the terminal.


### Option 1: Manual Connection (Temporary)
You must run this combined command every time you open a new terminal window.
``` bash
ssh -T git@github.com
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_tiger
ssh -T git@github.com
```

### Option 2: Config Connection (Permanent)
Set this up once, and it works automatically forever.
1. Open and create the config file:
``` bash
vim ~/.ssh/config # or create it manually under ~/.ssh folder
```
2. Paste the following settings:
```
host github.com
    HostName github.com
    User git
    IdentityFile ~/.ssh/id_tiger
    AddKeysToAgent yes
```

3. Test the connection:
``` bash
ssh -T git@github.com
```