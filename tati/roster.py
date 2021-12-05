import datetime as dt

# 1) put roster.jpeg through https://www.onlineocr.net/


base = dt.datetime(2021, 12, 6)

raw = """Mon 14-21:15
Tue
Wed 14-21:15
Thu
Fri 15-21:15
Sat 9:30-17:45
Sun 6:45-15
"""


def parse_hour(text):
    if ":" in text:
        hour = int(text.split(":")[0])
        minute = int(text.split(":")[-1])
    else:
        hour = int(text)
        minute = 0

    return hour, minute


start = []
end = []

for idx, line in enumerate(raw.split("\n")[:-1]):
    if len(line.split(" ")) < 2:
        continue

    day = base + dt.timedelta(days=idx)
    t0 = line.split(" ")[-1].split("-")[0]
    t1 = line.split(" ")[-1].split("-")[-1]

    hour, minute = parse_hour(t0)
    start.append(dt.datetime(day.year, day.month, day.day, hour, minute, 0))

    hour, minute = parse_hour(t1)
    end.append(dt.datetime(day.year, day.month, day.day, hour, minute, 0))
