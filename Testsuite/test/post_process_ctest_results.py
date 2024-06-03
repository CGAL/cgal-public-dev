import sys
import io
import re
import os

input_report_file=sys.argv[1]
report_file_name=sys.argv[2]
global_report_file_name=sys.argv[3]
config_regex=re.compile('(.*Configuring (examples|demo|test)*( in )*(test\/|examples\/|demo\/)*)((?!done)\w+)')
demo_regex=re.compile('.*in demo\/')
examples_regex=re.compile('.*in examples\/')
Separator = "------------------------------------------------------------------"


#open the Installation report
#For each NAME, check if NAME is a directory. If not, create one, create a 
#text report, and write everything that is in the report until the next NAME
#in it. Then, add 'NAME r' in the global report. This should allow to get all 
#the NOTICE and other info explaining why the configuration is skipped.

name=""
is_writing=False
is_ignored=False
position = 0
lines_to_write = []
installation_cmake_logs = []



def find_third_separator(contents):
    separator_count = 0
    for i, line in enumerate(contents):
        if line.strip() == Separator:
            separator_count += 1
            if separator_count == 3:
                return i
    return len(contents) + 2

def find_last_separator(contents):
    for i, line in enumerate(contents):
        if line.strip() == Separator:
            position = i
    return position

with open ("{dir}/{file}".format(dir="Installation",file=report_file_name), "r") as file:
  contents = file.readlines()
position = find_last_separator(contents)
for i, line in enumerate(contents):
    if i > position:
        installation_cmake_logs.append(line)
    if line.strip() == "== Generating build files for tests ==":
        break
contents = []

global_report = open(global_report_file_name, "a+")
with open(input_report_file, "rt") as test_report:
    for myline in test_report:
        match = config_regex.match(myline)
        if is_writing:
            if match:
                is_writing = False
                if lines_to_write:
                    file_path = "{dir}/{file}".format(dir=name, file=report_file_name)
                    if os.path.exists(file_path):
                        with open(file_path, "r") as file:
                            contents = file.readlines()
                    else:
                        contents = []

                    position = find_third_separator(contents)

                    if Separator + "\n- CMake Results \n" + Separator not in lines_to_write:
                        lines_to_write.insert(0,Separator + "\n- CMake Results \n" + Separator + "\n")
                    if Separator + "\n- CMake Logs \n" + Separator not in contents:
                        contents.insert(position - 1, Separator + "\n- CMake Logs \n" + Separator + "\n\n")
                    for log in installation_cmake_logs:
                        contents.insert(position, log)
                        position += 1
                    lines_to_write.insert(0, "\n")
                    contents[position:position] = lines_to_write

                    with open(file_path, "w") as file:
                        file.write("".join(contents))

                    lines_to_write = []

                if is_ignored:
                    is_ignored = False
            else:
                if myline.strip() != "":
                    lines_to_write.append(myline)
        if not is_writing:
          if match:
            name=match.group(0).replace(match.group(1), "")
            if demo_regex.match(myline):
              name="{str}_Demo".format(str=name)
            elif examples_regex.match(myline):
              name="{str}_Examples".format(str=name)
            elif name == "libCGAL":
              name="libCGAL_shared"
            elif name == "libCGAL_Core":
              name="libCGALCore_shared"
            elif name == "libCGAL_ImageIO":
              name="libCGALimageIO_shared"
            elif name == "libCGAL_Qt6":
              name="libCGALQt6_shared"
            if name=="incomplete":
              is_writing=False
              is_ignored=False
              continue
            else:
                if not os.path.isdir(name):
                    is_ignored = True
                    os.mkdir(name)
                    with open("{dir}/{file}".format(dir=name, file=report_file_name), "w") as test_report:
                        test_report.write(open("{}/../../../../../.scm-branch".format(os.getcwd()), 'r').read())
                else:
                    is_ignored = False
                is_writing = True

if is_writing:
    is_writing=False
    test_report.close()
    if is_ignored:
      print("{label} {result}".format(label=name, result='r'), file=global_report)
      is_ignored=False
global_report.close()
