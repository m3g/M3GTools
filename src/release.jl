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

run(`git add -A`)
run(`git commit -m "updated version file to $version"`)
run(`git tag -a $version -m "Release $version"`)
run(`git push origin master tag $version`)

range = "$(tags[length(tags)])...$version"
tagdiff = split(read(`git log --pretty=oneline $range`,String))
println(tagdiff)

#println("----------------------")
#println("CHANGE LOG:")
#println("----------------------")
#println(tagdiff)
#for line in eachline(tagdiff)
#  println(line)
#end
#println("----------------------")


