#!/bin/bash

# This proxy script is designed to be symlinked into
# a .git/hooks directory as one or more of the regular
# git hook names. It looks for scripts of the format
# <hookname>-<whatever> (e.g. post-receive-01-foobar)
# and invokes them with whatever arguments and stdin
# were originally passed to this script. This allows
# multiple different hook scripts to all be installed
# for the same Git hook event.
#
# Exit statuses are aggregated such that if any script
# exits nonzero, this proxy will also exit nonzero.
#
# ----------------------------------------------------
#
# Example usage:
#
#  $ ln -s githook-proxy.sh <somerepo>.git/hooks/post-receive
#  $ touch <somerepo>.git/hooks/post-receive-01-do-something-fun
#  $ touch <somerepo>.git/hooks/post-receive-02-do-other-stuff
#  $ chmod +x <somerepo>.git/hooks/post-receive-*
#
# Both post-receive-01-do-something-fun and post-receive-02-do-other-stuff
# would then be called with the same stdin and arguments as if
# they were a directly installed hook named 'post-receive'. If either
# hook exited nonzero, the proxy would exit nonzero.
#
# ----------------------------------------------------

set -u
set -e

stdin=$(cat)
me=$(readlink -e "$0") # Actual script path
hooktype=$(basename "$0") # what hook was invoked
hookdir=$(dirname "$0") # where we should look for hooks

# Default to exit 0 (e.g. if no hooks are run)
# If any hook exits nonzero, we'll exit nonzero
exitstatus=0

while read HOOK; do

	# Make sure we don't somehow recurse
	HOOK="$(readlink -e "$HOOK")"
	if [ "$HOOK" == "$me" ]; then continue; fi

	# Only try to run executable hooks
	if [ ! -x "$HOOK" ]; then continue; fi

	echo "Running $(basename "$HOOK") ..."

	# Invoke each found hook with identical
	# standard input and arguments as we were
	# called with.
	if ! echo "$stdin" | "$HOOK" "$@"; then
		exitstatus=1
	fi

done < <(find "$hookdir" -maxdepth 1 -name "$hooktype-*" | sort)

exit "$exitstatus"
