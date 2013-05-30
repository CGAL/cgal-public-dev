#!/bin/bash
PWD=$(pwd)
/var/lib/gems/1.8/bin/git-commit-notifier hooks/git-commit-notifier/${PWD:9}-grouped.yml
