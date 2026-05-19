from jarvis.db.jsonutils import loadjson
import glob

for i in glob.glob("*.json"):
    d = loadjson(i)
    fname = i.split(".json")[0] + ".csv"
    f = open(fname, "w")
    f.write("jid,formula,max_voltage,voltages\n")
    for j in d:
        line = (
            str(j["jid"])
            + ","
            + str(j["formula"])
            + ","
            + str(j["max_voltages"])
            + ","
            + str(j["voltages"])
            + "\n"
        )
        f.write(line)
    f.close()
