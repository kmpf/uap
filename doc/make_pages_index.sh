#! /bin/bash -

cd "$1"
touch index.html

printf '
<html>
<body>
<h2>UAP Sphinx Docu per branch</h2>
<p>
' >> index.html

for file in $(ls -1); do
    [[ "$file" == "index.html" ]] && continue
    printf '<li><a href=%s>%s</a></li>' "$file" "$file" >> index.html
done

printf '
</p>
</body>
</html>
' >> index.html
