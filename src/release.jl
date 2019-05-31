
using Dates
tags = split(read(`git tag`,String))
version = "$(Dates.year(Dates.today())-2000).$(Dates.dayofyear(Dates.today()))"
i = 1
while version in tags 
  global i, version
  i = i + 1
  version = version*".$i"
end
println("Will create version: $version")

run(`git add -A`)
run(`git commit -m "updated version file to $version"`)
run(`git tag -a $version -m "Release $version"`)
#run(`git push origin master tag $version`)

range = read(`git tag '|' tail -n 2 '|' xargs '|' sed 's! !...!'`)
tagdiff = read(`git log --pretty=oneline $range '|' awk '{$1=""; print "-"$0}'`)

println("----------------------")
println("CHANGE LOG:")
println("----------------------")
for line in eachline(tagdiff)
  println(line)
end
println("----------------------")


