#!/usr/bin/awk -f
awk 'NR==1 { next }
       { T=0
          for(N=2; N<=NF; N++) T+=$N;
          T/=(NF-1)
          print $1, T }' $1
