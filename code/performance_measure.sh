#!/bin/bash

if [ "$#" -lt 6 ]; then
  echo "Usage: $0 start end step epsilon name algoflag [-a angle] [-l x y] [-d dimension] [-n batchsize]"
  exit 1
fi

start=$1
end=$2
step=$3
epsilon=$4
name=$5
algoflag=$6

angle="180"
length_min="1"
length_max="50"
dimension="2"
batchsize="20"

while [[ "$#" -gt 6 ]]; do
  case $7 in
  -a)
    if [[ "$8" =~ ^[0-9]+$ ]] && [ "$8" -ge 0 ] && [ "$8" -le 180 ]; then
      angle=$8
      shift 2
    else
      echo "Error: -a must be followed by a number between 0 and 180."
      exit 1
    fi
    ;;
  -l)
    if [[ "$8" =~ ^[0-9]+$ ]] && [[ "$9" =~ ^[0-9]+$ ]]; then
      length_min=$8
      length_max=$9
      shift 3
    else
      echo "Error: -l must be followed by two numbers."
      exit 1
    fi
    ;;
  -d)
    if [[ "$8" =~ ^[0-9]+$ ]] && [ "$8" -ge 1]; then
      dimension=$8
      shift 2
    else
      echo "Error: -d must be followed by one number >= 1"
      exit 1
    fi
    ;;
  -n)
    if [[ "$8" =~ ^[0-9]+$ ]] && [ "$8" -ge 1]; then
      batchsize=$8
      shift 2
    else
      echo "Error: -n must be followed by one number >= 1"
      exit 1
    fi
    ;;
  *)
    echo "Unknown option: $7"
    exit 1
    ;;
  esac
done

main_dir="data/$name"
if [ -d "$main_dir" ]; then
  echo "Error: Directory $main_dir already exists."
  exit 1
fi

if [ ! -d "build" ]; then
  echo "Error: build directory does not exist"
  exit 1
fi

if [ ! -f "build/polyline" ]; then
  echo "Error: project has not been built yet"
  exit 1
fi

mkdir -p "$main_dir"

for ((i = start; i <= end; i += step)); do
  sub_dir="$main_dir/$i"
  mkdir -p "$sub_dir"

  datagen_cmd="./build/datagen $i $dimension -f $sub_dir/p -n $batchsize -a $angle -l $length_min $length_max"
  echo "Generate polyline test cases for length $i polylines"
  $datagen_cmd
done

mkdir -p "measurements/$name"
for ((i = start; i <= end; i += step)); do
  sub_dir="$main_dir/$i"

  polyline_cmd="./build/polyline $sub_dir $epsilon $algoflag"
  echo "Perform algorithm on $sub_dir"
  $polyline_cmd
  mv "measurements/$i/"* "measurements/$name/"
  rmdir "measurements/$i"
done

if [ -d "visualization" ]; then
  rm -r visualization
fi

echo -e "{\n  \"dimension\": $dimension,\n  \"max_angle\": $angle,\n  \"min_length\": $length_min,\n  \"max_length\": $length_max,\n  \"epsilon\": $epsilon,\n  \"algorithm_flag\": \"$algoflag\",\n  \"batchsize\": $batchsize\n}" >"measurements/$name/specification.json"
