#!/bin/bash
PWD=$(pwd)
/var/local/git-commit-notifier/bin/git-commit-notifier hooks/git-commit-notifier/${PWD:9}.yml
