#!/bin/bash

# Copyright (c) 2023 Stephen Wasilewski, EPFL
# =======================================================================
# This program is free software: you can redistribute it and/or
# modify it under the terms of theGNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.
# =======================================================================

help()
{
    echo "Usage: extract_data.sh [ -v LEFT BOTTOM WIDTH HEIGHT]
                       [ -h LEFT BOTTOM WIDTH HEIGHT]
                       [ -r LEFT BOTTOM WIDTH HEIGHT]
                       [ -c X Y Radius]
                       image_file
    -v: vertical section (average horizontally, sorted top -> bottom)
    -h: horizontal section (average vertically, sorted left -> right)
    -r: rectangular region (sorted highest -> lowest)
    -c: circular region, specify by pixel center and radius (in degrees)
        (sorted highest -> lowest)
    -i: output cropped hdr instead of data
    -rgb: output rgb instead of luminance
    
    only performs one operation, uses last options specified. Images must be square and have a valid VIEW= in header
    otherwise solid angle will revert to pixel
    prints cumulative solid angle, then pixel value, then cumulative average, .
    either omega L L_a or omega L R G B L_a R_a G_a B_a if -rgb.
"
    exit 2
}

while :
do
  case "$1" in
    -i )
      data="i"
      shift 1;
      ;;
    -rgb )
      rgb="rgb"
      shift 1;
      ;;
    -v )
      A=($2 $3 $4 $5)
      type="v"
      shift 5
      if [[ ${#A[@]} -ne 4 ]]; then
          echo "*** bad number of options for -v ***"
          echo
          help
      fi
      ;;
    -h )
      A=($2 $3 $4 $5)
      type="h"
      shift 5
      if [[ ${#A[@]} -ne 4 ]]; then
          echo "*** bad number of options for -h ***"
          echo
          help
      fi
      ;;
    -r )
      A=($2 $3 $4 $5)
      type="r"
      shift 5
      if [[ ${#A[@]} -ne 4 ]]; then
          echo "*** bad number of options for -r ***"
          echo
          help
      fi
      ;;
    -c )
      A=($2 $3 $4)
      type="c"
      shift 4
      if [[ ${#A[@]} -ne 3 ]]; then
          echo "*** bad number of options for -c ***"
          echo
          help
      fi
      ;;
    -* )
      help
      ;;
    *)
      img=$1
      break;
      ;;
  esac
done

if [ -z "$img" ]; then
    echo "*** missing image file (or bad # of options) ***"
    echo
    help
fi

basename="${img%.*}"
pr=($(getinfo -d "$img"))
xres=${pr[4]}

amg="$basename"_ang.hdr

hasview=$(getinfo $img | grep VIEW)
if [ -z "$hasview" ]
then
   pcomb -e 'lo=1/(xres*yres)' $img > "$amg"
else
  pcomb -e 'lo=S(1)' $img > "$amg"
fi


if [ "$rgb" == "rgb" ]; then
    pval="-o -h -H -d"
    
    if [ "$type" == "c" ]; then
        rcal='grey(r,g,b):.265074126*r+.670114631*g+.064811243*b;cond=$1;$1=grey($1,$2,$3)*179;$2=$1*179;$3=$2*179;$4=$3*179;$5=$4'
    else
        rcal='grey(r,g,b):.265074126*r+.670114631*g+.064811243*b;$1=grey($1,$2,$3)*179;$2=$1*179;$3=$2*179;$4=$3*179;'
    fi
    rcalf='$1=$5*recno;$2=$6;$3=$7;$4=$8;$5=$9;$6=$1;$7=$2;$8=$3;$9=$4'
else
    pval="-o -b -h -H -d"
    if [ "$type" == "c" ]; then
        rcal='cond=$1;$1=$1*179;$2=$2'
    else
        rcal='$1=$1*179;'
    fi

    rcalf='$1=$2*recno;$2=$3;$3=$1'
fi

if [ "$type" == "v" ]; then
    if [ "$data" == "i" ]; then
        pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]}
        exit 0
    fi
    datac=$((${A[2]}*${A[3]}))
    pcompos -x ${A[2]} -y ${A[3]} $amg -${A[0]} -${A[1]} | pvalue -o -b -h -H -d | total -${A[2]} -m  > "$basename"_ang.txt
    pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]} | pvalue $pval | total -${A[2]} -m  | rcalc -e ''$rcal'' | rlam - "$basename"_ang.txt > "$basename"_tmp.txt
elif [ "$type" == "h" ]; then
    if [ "$data" == "i" ]; then
        pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]}
        exit 0
    fi
    datac=$((${A[2]}*${A[3]}))
    pcompos -x ${A[2]} -y ${A[3]} $amg -${A[0]} -${A[1]} > "$basename"_ang2.hdr
    pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]} > "$basename"_tmp.hdr
    protate "$basename"_ang2.hdr | pvalue -o -b -h -H -d | total -${A[3]} -m  > "$basename"_ang.txt
    protate "$basename"_tmp.hdr | pvalue $pval | total -${A[3]} -m  | rcalc -e ''$rcal'' | rlam - "$basename"_ang.txt > "$basename"_tmp.txt
    rm "$basename"_ang2.hdr "$basename"_tmp.hdr
elif [ "$type" == "r" ]; then
    if [ "$data" == "i" ]; then
        pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]}
        exit 0
    fi
    pcompos -x ${A[2]} -y ${A[3]} $amg -${A[0]} -${A[1]} | pvalue -o -b -h -H -d > "$basename"_ang.txt
    pcompos -x ${A[2]} -y ${A[3]} $img -${A[0]} -${A[1]} | pvalue $pval | rcalc -e ''$rcal'' | rlam - "$basename"_ang.txt | sort -g -r  > "$basename"_tmp.txt
elif [ "$type" == "c" ]; then
    pr=($(rcalc -n -e 'ppd='"$xres"'/'"${vv[1]}"';prad=ppd*'"${A[2]}"';iframe=ceil(prad+1)*2;$1=ppd;$2=prad;iprad=ceil(prad);$3=iprad;$4=iframe;$5='"${A[0]}"'-iprad;$6='"${A[1]}"'-iprad'))
    if [ "$data" == "i" ]; then
        pcompos -x ${pr[3]} -y ${pr[3]} "$img" -${pr[4]} -${pr[5]} | pcomb -e 'r='"${pr[1]}"';ir='"${pr[2]}"';c(i)=if((x-ir)^2+(y-ir)^2-r^2,0,i);ro=c(ri(1));go=c(gi(1));bo=c(bi(1))' -
        exit 0
    fi
    pcompos -x ${pr[3]} -y ${pr[3]} "$amg" -${pr[4]} -${pr[5]} | pvalue -o -b -h -H -d > "$basename"_ang.txt
    pcompos -x ${pr[3]} -y ${pr[3]} "$img" -${pr[4]} -${pr[5]} | pcomb -e 'r='"${pr[1]}"';ir='"${pr[2]}"';c(i)=if((x-ir)^2+(y-ir)^2-r^2,0,i);ro=c(ri(1));go=c(gi(1));bo=c(bi(1))' - | pvalue $pval | rlam - "$basename"_ang.txt | rcalc -e ''$rcal'' | sort -g -r  > "$basename"_tmp.txt
fi

total -r -m -1 "$basename"_tmp.txt | rlam - "$basename"_tmp.txt | rcalc -e ''$rcalf''
 rm "$basename"_tmp.txt "$basename"_ang.txt "$amg"
