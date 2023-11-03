(k:recurse #P"include/" #P"src/")

(k:systems :cando-user)

(k:sources :install-bin
           #~"cando-user-install"
           #~"cando-apptainer-install"
           #~"scando-apptainer-install")
