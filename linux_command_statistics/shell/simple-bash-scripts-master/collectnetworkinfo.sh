#!/bin/bash

IP4FW=/sbin/iptables
IP6FW=/sbin/ip6tables
LSPCI=/usr/bin/lspci
ROUTE=/sbin/route
NETSTAT=/bin/netstat
LSB=/usr/bin/lsb_release

## files ##
DNSCLIENT="/etc/resolv.conf"
DRVCONF="/etc/modprobe.conf"
NETALIASCFC="/etc/sysconfig/network-scripts/ifcfg-eth?-range?"
NETCFC="/etc/sysconfig/network-scripts/ifcfg-eth?"
NETSTATICROUTECFC="/etc/sysconfig/network-scripts/route-eth?"
SYSCTL="/etc/sysctl.conf"

## Output file ##
OUTPUT="network.$(date +'%d-%m-%y').info.txt"

## Email info to?? ##
SUPPORT_ID="your_name@service_provider.com"

chk_root() {
	local meid=$(id -u)
	if [ $meid -ne 0 ]; then
		echo "You must be root user to run this tool"
		exit 999
	fi
}

write_header() {
	echo "---------------------------------------------------" >>$OUTPUT
	echo "$@" >>$OUTPUT
	echo "---------------------------------------------------" >>$OUTPUT
}

dump_info() {
	echo "* Hostname: $(hostname)" >$OUTPUT
	echo "* Run date and time: $(date)" >>$OUTPUT

	write_header "Linux Distro"
	echo "Linux kernel: $(uname -mrs)" >>$OUTPUT
	$LSB -a >>$OUTPUT

	[ -x ${HWINF} ] && write_header "${HWINF}"
	[ -x ${HWINF} ] && ${HWINF} >>$OUTPUT

	[ -x ${HWINF} ] && write_header "${HWINF}"
	[ -x ${HWINF} ] && ${HWINF} >>$OUTPUT

	write_header "PCI Devices"
	${LSPCI} -v >>$OUTPUT

	write_header "$IFCFG Output"
	$IFCFG >>$OUTPUT

	write_header "Kernel Routing Table"
	$ROUTE -n >>$OUTPUT

	write_header "Network Card Drivers Configuration $DRVCONF"
	[ -f $DRVCONF ] && grep eth $DRVCONF >>$OUTPUT || echo "Error $DRVCONF file not found." >>$OUTPUT

	write_header "DNS Client $DNSCLIENT Configuration"
	[ -f $DNSCLIENT ] && cat $DNSCLIENT >>$OUTPUT || echo "Error $DNSCLIENT file not found." >>$OUTPUT

	write_header "Network Configuration File"
	for f in $NETCFC; do
		if [ -f $f ]; then
			echo "** $f **" >>$OUTPUT
			cat $f >>$OUTPUT
		else
			echo "Error $f not found." >>$OUTPUT
		fi
	done

	write_header "Network Aliase File"
	for f in $NETALIASCFC; do
		if [ -f $f ]; then
			echo "** $f **" >>$OUTPUT
			cat $f >>$OUTPUT
		else
			echo "Error $f not found." >>$OUTPUT
		fi
	done

	write_header "Network Static Routing Configuration"
	for f in $NETSTATICROUTECFC; do
		if [ -f $f ]; then
			echo "** $f **" >>$OUTPUT
			cat $f >>$OUTPUT
		else
			echo "Error $f not found." >>$OUTPUT
		fi
	done

	write_header "IP4 Firewall Configuration"
	$IP4FW -L -n >>$OUTPUT

	write_header "IP6 Firewall Configuration"
	$IP6FW -L -n >>$OUTPUT

	write_header "Network Stats"
	$NETSTAT -s >>$OUTPUT

	write_header "Network Tweaks via $SYSCTL"
	[ -f $SYSCTL ] && cat $SYSCTL >>$OUTPUT || echo "Error $SYSCTL not found." >>$OUTPUT

	echo "The Network Configuration Info Written To $OUTPUT. Please email this file to $SUPPORT_ID."
}

chk_root
dump_info
