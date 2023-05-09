using Gumbo
using Cascadia
using Glob
using Franklin

function hfun_bar(vname)
  val = Meta.parse(vname[1])
  return round(sqrt(val), digits=2)
end

function hfun_m1fill(vname)
  var = vname[1]
  return pagevar("index", var)
end

function lx_baz(com, _)
  # keep this first line
  brace_content = Franklin.content(com.braces[1]) # input string
  # do whatever you want here
  return uppercase(brace_content)
end


get_html_files() = vcat(glob("*.html", "./_assets/notebooks"), glob("*.html", "./_assets/prerendered_notebooks"))

nb_name(rpath) = begin
  file = split(rpath, '/')[end]
  split(file,'.')[1]
end

function parse_tags_from_html(filepath::String)
  html_content = read(filepath, String)
  html = parsehtml(html_content)
  tags_selector = Selector("meta[property='og:article:tag']")
  tags_nodes = eachmatch(tags_selector, html.root)
  return map(x -> get(x.attributes, "content", ""), tags_nodes)
end

function parse_title_from_html(filepath::String)
  html_content = read(filepath, String)
  html = parsehtml(html_content)
  # Select the title element using Cascadia
  title_selector = Selector("head title")
  title_element = eachmatch(title_selector, html.root)

  # Extract the title text
  if !isempty(title_element)
    title = nodeText(title_element[1])
  else
    title = ""
  end
  return title
end



function hfun_taglist_pluto()::String
  tag = locvar(:fd_tag)::String

  c = IOBuffer()
  write(c, "<ul>")

  rpaths = globvar("fd_tag_pages")[tag]
  sorter(p) = begin
    pvd = pagevar(p, "date")
    if isnothing(pvd)
      return Date(Dates.unix2datetime(stat(p).ctime))
    end
    return pvd
  end
  sort!(rpaths, by=sorter, rev=true)



  for rpath in rpaths
    title = parse_title_from_html(replace(rpath, "assets" => "_assets"))

    if isempty(title)
      title = splitext(splitpath(rpath)[end])[1]
      # title = join(nb_name(rpath),".jl")
    end

    url = get_url(rpath)

    write(c, "<li><a href=\"/$rpath\">$title</a></li>")
  end
  write(c, "</ul>")
  return String(take!(c))
end




"""
Working!
Add tags to Franklin variable.
"""
function hfun_get_pluto_tags()
  html_files = get_html_files()

  # Grab all tag names first
  tags = []
  rpaths = []
  for html_file in html_files
    push!(rpaths,replace(html_file, "_assets" => "assets"))
    itag = Franklin.refstring.(parse_tags_from_html(html_file))
    push!(tags, itag)
  end
  frmt_tags = Franklin.DTAG(zip(rpaths,map(Set,tags)))
  Franklin.set_var!(Franklin.GLOBAL_VARS, "fd_page_tags", frmt_tags; check=false)

  return ""
end


"""
Use to add tags on side of page.
"""
function hfun_get_tags()

  tags = Set{String}()

  html_files = get_html_files()

  for html_file in html_files
    file_tags = parse_tags_from_html(html_file)
    union!(tags, file_tags)
  end

  # Sort and collect the tags
  sorted_tags = sort(collect(tags))


  # Generate the HTML for the tag list
  html_buffer = IOBuffer()
  for tag in sorted_tags
    frmt_tag = join([titlecase(word) for word in split(tag)], " ")
    write(
      html_buffer,
      """
<li class="tag-item">
  <a class="tag-link" href="/tag/$(Franklin.refstring(tag))/">$(frmt_tag)</a>
</li>
"""
    )
  end

  return String(take!(html_buffer))
end
