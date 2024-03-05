import random
import matplotlib.pyplot as plt
import numpy as np

buffer_size = 60
#arrival_rate = 1.0 / 5.1e-9
#server_rate = 1.0 / 69.8e-9

arrival_rate = 1.0
server_rate = 0.1

#messages = 2000
messages = 200

# Model arrivals as Poisson, with exponentially distributed interarrival and
#  service times
interarrival_times = [random.expovariate(arrival_rate) for i in range(0, messages)]
service_times = [random.expovariate(server_rate) for i in range(0, messages)]

#plt.figure()
#plt.hist(interarrival_times, bins=int(np.floor(np.sqrt(messages))))

#plt.figure()
#plt.hist(service_times, bins=int(np.floor(np.sqrt(messages))))

# Model arrivals as bursty, with quick bursts of messages followed by longer
#  gaps

interarrival_times = [5.1e-9] * 6  #ns
interarrival_times.append(random.uniform(30e-9, 90e-9))  # ns
#interarrival_times.append(300e-9)  # ns
interarrival_times = interarrival_times * 100
interarrival_times = interarrival_times[0:messages]
print(interarrival_times)

# Model service as uniformly distributed
#service_times = [random.uniform(2.31E-8, 1.31E-7) for _ in range(0, messages)]
service_times = [70e-9] * messages

# ** End of arrival and service models **


times = [0,]
queue_sizes = [0,]

m = 0
queue_len = 0

t = 0
updates = []
for i in interarrival_times:
    assert(i > 0)
    t += i
    updates.append((t, +1),)
updates = sorted(updates, key=lambda u: u[0], reverse=True)


event_count = 0
while len(updates) > 0:
    event_count += 1
    if event_count % 100 == 0:
        print(event_count)
    u = updates.pop()
    t = u[0]
    if u[1] == +1:  # insert
        if queue_len < buffer_size:
            queue_len += 1
            assert(len(service_times) > 0)
            if queue_len == 1:  # If this is the head of the queue, schedule service
                updates.append((t + service_times.pop(), -1))

    elif u[1] == -1:  # service
        queue_len -= 1
        if queue_len > 0:
            updates.append((t + service_times.pop(), -1))

    times.append(t)
    queue_sizes.append(queue_len)
    updates = sorted(updates, key=lambda u: u[0], reverse=True)


plt.figure()
plt.plot(times, queue_sizes, '-')

plt.show()
