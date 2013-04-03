#!/usr/bin/env ruby1.9.1
require 'rubygems'
require 'RedCloth'
require 'json'
require 'digest/sha1'
require 'fileutils'
require 'base64'

LEGAL_FUNCTIONS = ['includeJsonSnippet', 'includeSnippet', 'toc', 'dot']

FileUtils::mkdir('cache') if !File::directory?('cache')

def includeJsonSnippet(path)
    return "<notextile><listing>" + JSON.pretty_generate(JSON.parse(File::read(path))) + "</listing></notextile>"
end

def includeSnippet(path)
    return "<notextile><listing>" + File::read(path) + "</listing></notextile>"
end

def toc(toc_levels)
    result = ''
    count = 1
    counter = []
    $contents.scan(/^h(\d+)\.\s+(.+)$/).each do |hit_x|
        level = hit_x[0].to_i
        caption = hit_x[1].to_s

        while counter.size < level
            counter << 0
        end

        while counter.size > level
            counter.pop()
        end

        counter[counter.size - 1] = counter.last + 1

        enum = ''
        if (toc_levels.include?(level))
            if level < 4
                enum = counter[toc_levels.first - 1, counter.size].join('.') + '&nbsp;&nbsp;'
            end
#             result += "#{'*' * (level - toc_levels.first + 1)} \"#{enum}#{caption}\":#h#{count}\n"
            tocCaption = caption.dup
            if tocCaption.include?('~')
                tocCaption = tocCaption[0, tocCaption.index('~')]
            end
            result += "#{'*' * (level - toc_levels.first + 1)}#{level == 1 ? '(first)': '(later)'} #{enum}<a href='#h#{count}'>#{tocCaption}</a>\n"
        end
        $contents.sub!(/^h(\d+)\.\s+(.+)$/, "h#{hit_x[0].to_i}(#h#{count}). #{enum}#{caption.gsub('~', ' ')}")
        count += 1
    end
    return "<div class='toc'>\n#{result}\n</div>"
end

def dot(path)
    outPath = path + '.svg'
    system("dot -Tsvg -o \"#{outPath}\" \"#{path}\"")
    return "<img src='#{outPath}' />"
end

template = File::read('template.xhtml')
$contents = File::read('documentation.textile')

# scan for /* ... */ comments
while $contents.index('/*') != nil
    i0 = $contents.index('/*')
    i1 = $contents.index('*/')
    break unless i1
    $contents[i0..(i1 + 1)] = ''
end

# scan for <dot> ... </dot>
while $contents.index('<dot>') != nil
    i0 = $contents.index('<dot>')
    i1 = $contents.index('</dot>')
    break unless i1
    hit = $contents[i0 + 5, i1 - i0 - 5].strip
    sha1 = Digest::SHA1.hexdigest(hit)
    cachePath = "cache/#{sha1}.png"
    png = ''
    if File::exists?(cachePath)
        png = File::read(cachePath)
    else
        IO::popen("dot -Tpng -o /dev/stdout /dev/stdin", 'r+') do |io|
            io.puts hit
            io.close_write()
            png = io.read
        end
        File::open(cachePath, 'w') do |f|
            f.puts png
        end
    end
    $contents[i0..(i1 + 5)] = "<notextile><img src='data:image/png;base64,#{Base64.encode64(png)}' /></notextile>\n"
end


# scan for functions
regex = /#\{([^\}]+)\}/
$contents.scan(regex).each do |hit_x|
    hit = hit_x.first.to_s
    ok = false
    LEGAL_FUNCTIONS.each do |x|
        if hit.index(x) == 0
            ok = true
        end
    end
    if ok
        $contents.sub!(regex, eval(hit))
    else
        puts "Error: Illegal function call (#{hit})."
    end
end

$contents.gsub!("<listing>\n", '<listing>')
$contents.gsub!('<listing>', '<notextile><listing>')
$contents.gsub!('</listing>', '</listing></notextile>')

$contents = RedCloth.new($contents).to_html

template.sub!('#{TEXTILE}', $contents)
# remove all line breaks, we don't need them
template.gsub!('<br />', '')

File::open('documentation.xhtml', 'w') do |f|
    f.puts template
end
