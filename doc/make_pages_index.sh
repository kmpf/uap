#! /bin/bash -

cd "$1"

printf '
<html>
<head>
	<link rel="stylesheet" href="styles.css">
	<link rel="apple-touch-icon" sizes="180x180" href="/apple-touch-icon.png">
	<link rel="icon" type="image/png" sizes="32x32" href="/favicon-32x32.png">
	<link rel="icon" type="image/png" sizes="16x16" href="/favicon-16x16.png">
	<link rel="manifest" href="/site.webmanifest">
</head>
<body>
	<div>
	<h2>UAP Sphinx docu per branch</h2>
	<p>
' > index.html

for file in $(ls -1); do
    [[ "$file" == "index.html" ]] && continue
    [[ ! -d "$file" ]] && continue
    printf '		<li><a href=%s>%s</a></li>\n' \
		"$file" "$file" >> index.html
done

printf '
	</p>
	</div>
</body>
</html>
' >> index.html
