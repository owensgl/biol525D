---
title: "Key Authentication: configure your ssh session"
layout: page
---

This page contains instructions to connect rapidly to your account. The exercise is part of [Topic 1](./). It assumes you've completed the steps: [generate a key](./generate_a_key) and [configure ssh-agent](./configure_ssh_agent).

> *Note: watch out for placeholders* Replace the placeholders `serveruser`, `serverpass`, and `serverhost`, with your assigned username, password, and server address (IP).

Time to complete: 5 minutes

### Instructions for MacOS ###

1. Edit (or create) file `~/.ssh/ssh_config`, and add this information to it, replacing placeholders with your data, (and then save the file):

       # Allows login to your assigned server by simply doing: ssh biol525
       Host biol525
       HostName serverhost
       IdentityFile ~/.ssh/biol525D
       User serveruser
       ForwardAgent yes
       ForwardX11 yes

   Explanation:
   
     - `Host` setting is a friendly name that you choose. You will use it as a shorthand alias when connecting: i.e. `ssh biol525`.
	 - `HostName` is the address of your assigned server (e.g. 123.123.123.123).
     - `IdentityFile` is the path to the private key to use (it will use that one from your agent).
	 - `User` is your username _on the server_.
	 - `ForwardAgent` allows you to connect to other servers _from that server_ while relying on the local agent. This is useful if you hop on from one server to another using one set of keys. If you configure your public key in your github accout, for instance, you would be able to issue git commands without passwords, _on the server_. If you don't trust the administrators of the server, you would choose "no".
     - `ForwardX11` allows running graphical applications installed on the server, but viewing the window _from your computer_. Again, if you don't trust the admins, set it to "no".
	 - If you wish to know more about this file, see the manual page: `man ssh_config`

1. Configuring cyberduck.

   No config needed! Once you've configured your `ssh_config` like
   described above, then you should be able to connect to your server
   using the short alias you configured, e.g. `biol525`. Cyberduck will
   also automatically rely on your `ssh-agent` to provide the
   necessary credentials.

1. The remaining step is to configure your account on the server to recognize the public key that is loaded in your ssh-agent's keyring. This is done by logging in _once_ using your password, and then adding the public key to a recognized filename. There is a helper command that does it for you:


    # what this does:
	#   - logs in to the server as serveruser
	#   - appends the given public key to ~/.ssh/authorized_keys.
	# replace serveruser and serverhost with your assigned credentials
    ssh-copy-id -i ~/.ssh/biol525D.pub serveruser@serverhost

	# it will ask you for your account password --
	# once this is configured, you won't need your password again to connect over ssh.

    > *Note:* The first time you connect to a server over ssh, it will
	>         ask you to authorize the server's public keys. The
	>         server is also using a public key to identify itself to
	>         your computer. Ideally, you know ahead of time which
	>         public key to expect. If you're interested you can read
	>         up on "man in the middle attack". If you're connecting
	>         to a server which is on your internal network, it's less
	>         of a concern. You mark the key as trusted, and the ssh
	>         tools will remember the key for that IP.

### Instructions ###
