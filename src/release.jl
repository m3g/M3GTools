using Dates
tags = split(read(`git tag`,String))
version_base = "$(Dates.year(Dates.today())-2000).$(Dates.dayofyear(Dates.today()))"
i = 1
version = version_base
while version in tags 
  global i, version, version_base
  i = i + 1
  version = version_base*".$i"
end
println("Will create version: $version")

# Update version file:

version_func = "
function version()
  version = \"$version\"
  println(\"Version: \$version\")
end
"

file = open("./version.jl","w")
write(file,version_func)
close(file)
run(`git add -A`)
run(`git commit -m "updated version file to $version"`)
run(`git tag -a $version -m "Release $version"`)
run(`git push origin master tag $version`)

range = "$(tags[length(tags)])...$version"
tagdiff = split(read(`git log --pretty=oneline $range`,String),'\n')

println("----------------------")
println("CHANGE LOG:")
println("----------------------")
for line in tagdiff
  if length(line) < 1
    continue
  end
  num = split(line)[1]
  line = replace(line, num => "-")
  println(line)
end
println("----------------------")






