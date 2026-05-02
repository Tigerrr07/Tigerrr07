# Issue Summary Codex Extension Login Failure

When using VS Code Remote SSH with the Codex Extension to log in to ChatGPT, the browser redirects to
http://localhost:1455/auth/callback and returns ERR_CONNECTION_REFUSED.

**Cause**

The login process starts a callback server on the remote machine at port 1455. On the local machine, this port is forwarded to a dynamically assigned port, such as 58972. However, the browser still uses localhost:1455, which does not match the actual local port, resulting in a connection refusal.

**Solution**

Update the callback URL in the browser to use the dynamically assigned local port instead of 1455. For example, replace localhost:1455 with localhost:58972.

**Summary**

This issue is caused by a mismatch between the remote port and the dynamically assigned local port. The fix is to use the correct local port in the browser callback URL.