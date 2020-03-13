with open("unimod.xml", encoding="utf8") as f:
    mod_list = []
    while True:
        line = f.readline()
        if line == "": break
        elif "<umod:mod title=" in line:
            _start = line.find('title="') + len('title="')
            _end = line.find('"', _start)
            title = line[_start:_end]
            print(line)
            while not "record_id" in line:
                line = f.readline()
            _start = line.find('record_id="') + len('record_id="')
            _end = line.find('"', _start)
            id = line[_start:_end]
            mod_list.append((title, id))
with open("unimod.py", "w") as f:
    for title, id in mod_list:
        f.write("'%s',%s\n"%(title, id))