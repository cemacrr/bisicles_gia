#!/usr/bin/awk -f
{if ($9 =="dt") print $2 " " ($11)}
