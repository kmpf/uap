#!/usr/bin/env bash

echo "[make_doc.sh] Creating documentation as HTML"
make html && echo "[make_doc.sh] HTML has been successfully created."

while [[ "$pdf" != "y" ]] && [[ "$pdf" != "n" ]]; do
    read -p "[make_doc.sh] Create PDF of documentation? [y/n] " pdf
done

if [[ "$pdf" == "y" ]] ; then
    echo "[make_doc.sh] Creating documentation as PDF"
    make latexpdf && echo "[make_doc.sh] PDF has been successfully created."
fi