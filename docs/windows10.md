---
layout: "default"
title: Windows10
---

<h1><p style="text-align: center">Windows 10</p></h1>

-----
<br>

Windows 10 Ubuntu subsystem app setup
=================

Install Ubuntu 18 LTS from the Microsoft app store

Launch app

If needed Enable-WindowsOptionalFeature.  Follow instructions provided by link.  Run Powershell as administrator.

Launch Ubuntu app.  You will be prompted to make a username and password.  Don't forget them.

Your home directory will be: ~$ /home/<username>

Download Linux Anaconda Installer.  https://www.anaconda.com/download/#linux , Right click on "64-Bit(86) Installer" link and "Copy link address".

wget copied link

`~$ bash ./Anaconda3*` and follow prompts (installation will take a while, if it stalls, press <enter>)

Follow the Linux/Mac instructions at https://usda-vs.github.io/vSNP/setup.html start by cloning vSNP repository.

A note about the subsystem...

From within the Ubuntu subsystem terminal the directory structure is as a familiar Linux structure.  Directories under root are as expected with any Linux environment.  /home, /usr, etc.  The Windows directory structure is under /mnt.  With the c drive as a subdirectory under mnt.  Therefore your Windows 10 Desktop may have a similar path. ~$ /mnt/c/Users/<windows username>/OneDrive/Desktop
