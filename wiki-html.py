# given an XML Wiki file and a search pattern, produce an HTML page presenting each instance of the pattern in context, with interactive elements allowing associated text to be grabbed

import sys
import re
from re import finditer

# some example patterns and files
#pattern = "[Uu]nicell|[Uu]ni-cell|[Ss]ingle-cell"
#pattern = "[Ff]ree-living|[Ff]reeliving"
#pattern = "[Ex]tremophil"
#pattern = "palm"
#pattern = "ulticell|ulti-cell"
#pattern = "arasit"
#pattern = "[Ss]essil"
#pattern = "[Ff]lagell"
#wikifile = "mt-wiki.xml"
#outfilename = "form-mt-flag.html"

pattern = sys.argv[1]        # pattern to find
wikifile = sys.argv[2]       # XML dump from Wikipedia
outfilename = sys.argv[3]    # name for HTML output

capturepattern = "("+pattern+")"
buffer = 200                 # how many characters to present either side of match

titletext = "Results for /"+pattern+"/ in "+wikifile
print(titletext)

# open Wikipedia XML
fp = open(wikifile, "r")

# open HTML file for output 
outfile = open(outfilename, "w")

# write HTML preamble
outfile.write("<html>\n<head>\n<title>"+titletext+"</title>\n")
outfile.write("<style type=\"text/css\">\nhtml * {font-family: Arial;}\ntable, th, td { border: 1px #8888FF solid; background-color: #EEEEFF; border-collapse: collapse; table-layout: fixed; word-wrap: break-word; }\n</style>\n")
outfile.write("</head><body><h1>"+titletext+"</h1>\n")
outfile.write("<table width='1100px'>\n")

counter = 0
# go through Wiki file lines
for line in fp.readlines():
    # grab page title if present
    if "<title>" in line:
        titleline = line
    # loop through matches to our search pattern on this line
    for match in finditer(pattern, line):
       # extract the text either side of the pattern
       index = match.span()[0]
       substr = line[(max(0,index-buffer)):(min(len(line),index+buffer))]
       pagename = titleline.split("title")[1][1:-2]
       # first HTML column: page title
       outfile.write("<tr><td>  "+pagename+"  </td>")
       # second HTML column: checkbox that populates the final text box with the page title
       outfile.write("<td><input type='checkbox' onclick=\"document.getElementById('id"+str(counter)+"').value = document.getElementById('id"+str(counter)+"').value + '"+pagename+"';\"> </td>")
       # replace markup for italics (species names) with quotes
       substr = substr.replace("]]","''")
       substr = substr.replace("[[","''")
       # get locations of quotes -- we'll put a check box after each one
       quotes = [m.start() for m in re.finditer("''", substr)]
       quotes.insert(0,0)
       quotes.insert(len(quotes),2*buffer)
       tmpstr = ""
       # loop through quote locations in string
       for i in range(1,len(quotes)):
           # pull the text since the last quote location and make all instances of the search pattern bold
           tmpsubstr = re.sub(capturepattern, r"<b>\1</b>", substr[(quotes[i-1]):quotes[i]])
           # add an HTML checkbox that will populate the final textbox with the text since the last quote
           tmpstr = tmpstr+tmpsubstr+"<input type='checkbox' onclick=\"document.getElementById('id"+str(counter)+"').value = document.getElementById('id"+str(counter)+"').value+'"+tmpsubstr.replace("'","").replace("<b>","").replace("</b>","")+",';\">"
       # third HTML column: text and checkboxes
       outfile.write("<td width = '600px'> "+tmpstr+" </td>")
       # fourth HTML column: final textbox
       outfile.write("<td> <input type='text' id='id"+str(counter)+"'> </td></tr>\n")
       counter=counter+1

# final HTML: button and Javascript to populate a big textarea with all the textbox contents
outfile.write("</table>")
outfile.write("<script>\nfunction summary() {\nvar x = document.getElementsByTagName('INPUT');\nvar y = [];\nfor (var cnt = 0; cnt < x.length; cnt++) {\n  if (x[cnt].type == 'text') y.push(x[cnt].value);\n}document.getElementById('summarybox').value = '"+titletext+"'+y;\n}\n</script>")
outfile.write("<br><input type='button' value='Produce summary' onclick=\"summary()\"> <br><textarea id='summarybox' rows='20' cols='100'></textarea>")
outfile.write("</body></html>")
outfile.close()



    
    
