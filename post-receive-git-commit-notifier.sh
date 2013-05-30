#!/bin/sh
PWD=$(pwd); echo ${PWD:9}
/var/lib/gems/1.8/bin/git-commit-notifier hooks/git-commit-notifier/${PWD:9}.yml
