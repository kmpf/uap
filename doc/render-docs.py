#!../python_env/bin/python
import textile

template = '''
<!DOCTYPE html>
<html lang='en'>
<head>
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>RNASeq pipeline</title>

<style type='text/css'>
body {
    font-family: sans;
    font-size: 90%;
}

code {
    background-color: #f8f8f8;
    border: 1px solid #aaa;
    border-radius: 4px;
    padding-left: 0.3em;
    padding-right: 0.3em;
    padding-top: 0.1em;
    padding-bottom: 0.1em;
    font-size: 115%;
}

pre {
    background-color: #f8f8f8;
    border: 1px solid #aaa;
    border-radius: 8px;
    padding: 0.5em;
    font-size: 115%;
}

pre code {
    background: 0;
    border: 0;
    border-radius: 0;
    padding: inherit;
    font-size: 100%;
}

p {
    line-height: 1.5em;
}

</style>
</head>

<body>
#{CONTENT}
</body>
</html>
'''

with open('documentation.textile', 'r') as f:
    doc = f.read()

with open('documentation.html', 'w') as f:
    html = template.strip().replace('#{CONTENT}', textile.textile(doc))
    f.write(html)