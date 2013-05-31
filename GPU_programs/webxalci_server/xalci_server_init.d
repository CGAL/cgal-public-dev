#!/bin/sh
# Start/stop the xalci daemon.
#

test -f /home/emeliyan/work/binaries/bin/xalci_server || exit 0

. /lib/lsb/init-functions

case "$1" in
start)	log_daemon_msg "Starting xalci server" "xalci"
	start-stop-daemon --quiet --start --startas /home/emeliyan/work/binaries/bin/xalci_server --name xalci_server_daemon

        log_end_msg $?
	;;
stop)	log_daemon_msg "Stopping xalci server" "xalci"
        start-stop-daemon --quiet --stop --startas /home/emeliyan/work/binaries/bin/xalci_server --name xalci_server_daemon
        log_end_msg $?
        ;;
restart) log_daemon_msg "Restarting xalci server" "xalci" 
        start-stop-daemon --quiet --stop --startas /home/emeliyan/work/binaries/bin/xalci_server --name xalci_server_daemon
        log_end_msg $?
        ;;
*)	log_action_msg "Usage: /etc/init.d/xalci_server {start|stop|restart}"
        exit 2
        ;;
esac
exit 0
