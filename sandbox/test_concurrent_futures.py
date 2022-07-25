import cdsapi
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings("ignore")
 
LEADTIMES = ["%d" % (l) for l in range(24, 1128, 24)]
YEARS = ["%d" % (y) for y in range(1999, 2019)]
 
 
def get_dates(start=[2019, 1, 1], end=[2019, 12, 31]):
    start, end = datetime(*start), datetime(*end)
    days = [start + timedelta(days=i) for i in range((end - start).days + 1)]
    dates = [
        list(map(str.lower, d.strftime("%B-%d").split("-")))
        for d in days
        if d.weekday() in [0, 3]
    ]
    return dates
 
 
DATES = get_dates()
 
def retrieve(client, request, date):
 
    month = date[0]
    day = date[1]
    print(f"requesting month: {month}, day: {day} /n")
    request.update({"hmonth": month, "hday": day})
    client.retrieve(
        "cems-glofas-reforecast", request, f"glofas_reforecast_{month}_{day}.grib"
    )
    return f"retrieved month: {month}, day: {day}"
 
 
def main(request):
    "concurrent request using 10 threads"
    client = cdsapi.Client()
    with ThreadPoolExecutor(max_workers=10) as executor:
        futures = [
            executor.submit(retrieve, client, request.copy(), date) for date in DATES
        ]
        for f in as_completed(futures):
            try:
                print(f.result())
            except:
                print("could not retrieve")
 
 
if __name__ == "__main__":
 
    request = {
        "system_version": "version_2_2",
        "variable": "river_discharge_in_the_last_24_hours",
        "format": "grib",
        "hydrological_model": "htessel_lisflood",
        "product_type": "control_reforecast",
        "hyear": YEARS,
        "hmonth": "",
        "hday": "",
        "leadtime_hour": LEADTIMES,
    }
 
    main(request)
