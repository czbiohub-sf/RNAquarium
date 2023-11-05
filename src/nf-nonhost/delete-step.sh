grep $1 reports/trace-*.txt | cut -f2 | xargs -I% sudo find work -path "work/%*" -print -delete
