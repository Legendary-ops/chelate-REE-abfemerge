#!/bin/bash

declare -a dirs=(
 "e7c0409463b96d94d80ca90df746039b"
 "52fc28fb6fbff59ce3344f165f21b2b7"
 "8441ba19781fa583a1efa684f5a3eb4c"
 "291310fe57deb69fed0c907b387713d2"
 "e00c009c603c1d0d34b7194600e2491e"
 "c2597a69b4c32058a05cdcda3aaabfb0"
 "2056cef137509b280074b13bfd6f9df5"
 "8f74e581b3f536457115771179ad01e7"
 "418de272c7313c1f780685002bf0a685"
 "73a8cb71c3ee810b8ff749edbaf905a3"
 "1b30f8cb02e1b4a03efe3a75e54d3564"
 "f6a9ddd1a3ff410abfa9181c2871a808"
 "863f080480527b20b6c984b8b6a41554"
 "a8fa8ffac0e09822c8348676c42bae6e"
 "f3f1c64e49072cd165f7a55f48179e07"
 "27277d2836a8b3219ecf93a8d232cffb"
 "9d7c7258fb94ecc17234b9ba628c83c7"
 "c91d8fde34cba56e18c291f0a131c34b"
 "b3711a7532fc86c587333d2da4a9ec6c"
)



destination_folder="000000_analysis/Nd"
shopt -s nullglob

for dir_name in "${dirs[@]}"; do
    if [ -d "$dir_name" ]; then # Check if the directory exists
        echo "Processing directory: $dir_name"
        files_to_copy="$dir_name/PRO_CANON_*.xvg"
        matching_files=( $files_to_copy )
        if [ ${#matching_files[@]} -gt 0 ]; then
            cp -v "${matching_files[@]}" "$destination_folder/"
        else
            echo "  No files matching 'PRO_CANON_*.xvg' found in $dir_name. Skipping copy for this directory."
        fi
    else
        echo "Directory not found: $dir_name. Skipping."
    fi
done

shopt -u nullglob
