import datetime, time

start_time = time.time()

sum = 0 
for i in range(10000000):
    sum += 1

end_time = time.time()
time_interval = end_time - start_time
cost_time_list = str(datetime.timedelta(seconds=time_interval)).split(".")
#print(result_list[0])
print(cost_time_list[0])