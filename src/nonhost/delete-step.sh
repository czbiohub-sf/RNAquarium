grep $1 reports/trace-*.txt | cut -f2 | xargs -I% find work -path "work/%*" -print -delete
