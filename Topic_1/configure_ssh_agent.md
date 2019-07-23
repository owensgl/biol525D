---
title: "Configuring the SSH Agent"
layout: page
---


This page contains instructions to load up a key into the ssh-agent. The exercise is part of [Topic 1](./). It assumes you have followed steps in [Generate a key](./generate_a_key).

Instructions: Read the general text, and then follow only the instructions that pertain to your operating system and software.

Time needed to complete: 5 minutes.


One way to think about the ssh-agent is that it is a program which has a keyring on which you add keys. When a program needs to log in to a server over SSH, it will check with the agent to see if it has a valid key for that server.

> ssh-agent is a program to hold private keys used for public key authentication (RSA, DSA). The idea is that ssh-agent is started in the beginning of an X-session or a login session, and all other windows or programs are started as clients to the ssh-agent program. Through use of environment variables the agent can be located and automatically used for authentication when logging in to other machines using ssh(1).
>
> (excerpt from the [manual page](https://linux.die.net/man/1/ssh-agent))

_Other Info_:
  - [wikipedia](https://en.wikipedia.org/wiki/Ssh-agent)



### Instructions for Mac or GNU/Linux ###

1. Open a terminal.

1. Check if your agent is already running:

        # lists identities currently loaded in the agent
        ssh-add -l

   - If it says "The agent has no identities." then it is running, but no keys are loaded in it yet. Skip step 3 (starting the agent).
   - If it lists one or more key fingerprints, such as

         4096 SHA256:W/q4wBqfovIxIL7iO6mnDM/IPl0pjW5SqbjwinYmu70 /home/foo/.ssh/id_rsa2 some name (RSA)

     then skip step 3 (starting the agent).

   - If it says It "cannot connect to the agent", or "Connection
     refused" then it needs to be started. This would be surprising,
     but not impossible. Typically the agent will be available from
     all terminals started locally in your desktop session.

1. Starting the agent: For whichever reason, _if it's not already running_, you can start an Agent with the following command:

        eval "$(ssh-agent -s)"

   *Note:* This will start the agent in the background of the _current_ shell session _only_. (Optional: To ensure your agent is running
           whenever you open a new terminal (or another tab), you may add the following to the end of file `~/.bashrc`):

    ```
	#
    # ssh-agent configuration
    #
    if [ -z "$(pgrep ssh-agent)" ]; then
        rm -rf /tmp/ssh-*
        eval "$(ssh-agent -s)" > /dev/null
    else
        export SSH_AGENT_PID=$(pgrep ssh-agent)
        export SSH_AUTH_SOCK=$(find /tmp/ssh-* -name agent.*)
    fi
	```

1. Add the key to your agent (you want the private key, not the .pub). You will be prompted for your passphrase. You will have to perform this operation once _per computer session_ (i.e. next time you reboot your computer).

       ssh-add ~/.ssh/biol525D


### Instructions for Windows: MobaXTerm ###

1. Open MobaXTerm settings. `Menu` -> `Settings` -> `Configuration`

1. Under the `SSH` tab, in the `SSH Agents` section, enable `Use Internal SSH Agent MobAgent`

1. In the list of keys to load on Startup, click `+` and add your `biol525D.ppk` file that you generated following the [previous instructions](./generate_a_key). You will be prompted to restart MobaXTerm.

1. When Moba restarts, it will prompt you for the passphrase for your key.

1. If you go back to that SSH configuration tab later, in the configuration, you should see your key in the list.


### Instructions for Windows: PuTTY ###

The PuTTY suite comes with Pagent.exe which runs in the background. When it is running, it will signal its presence with an icon in the system tray.

1. Run Pagent.exe. <kbd>Winkey</kbd>+<kbd>r</kbd>, type Pageant + <kbd>Enter</kbd>

1. Pageant runs silently, but it adds a little icon in your system tray -- a computer with a hat. Right click on the icon and choose "Add Key".

1. Select your `biol525D.ppk` file that you selected earlier, and enter its passphrase. You will have to perform this operation _once_ every time you login to your computer after boot. Optional: You may add Pageant to the list of programs to start automatically when you login to Windows.

1. If you don't remember if your key is loaded or not, you can do "List keys" in the right-click menu of Pageant.

### Next step ###

Follow instructions at [finalize tool config](./finalize_tool_config).
