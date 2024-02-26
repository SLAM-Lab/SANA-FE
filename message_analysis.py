import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import math
import os
import queue

np.set_printoptions(threshold=100000)
np.seterr(invalid='raise')

def hw_str_to_core(hw_str):
    tile, core = hw_str.split('.')
    core = int(core)
    tile = int(tile)
    return tile, core


def graph_index(x, y, link):
    assert(x < 8)
    assert(y < 4)
    assert(link < 12)
    return (x*4*12) + (y*12) + link


def reverse_graph_index(index):
    x = index // (4*12)
    y = (index % (4*12)) // 12
    link = index % 12
    assert(x < 8)
    assert(y < 4)
    assert(link < 12)
    return (x, y, link)


def calculate_a_k(k, arrival_rate, service_distribution_pdf):
    x = service_distribution_pdf[0]
    y = ((arrival_rate * x) ** k) / math.factorial(k)
    y = y * np.exp(-arrival_rate * x) * service_distribution_pdf[1]

    #print(f"arrival:{arrival_rate} k:{k} dist:{service_distribution_pdf[1]}")
    #print(f"x:{x} y:{y}")

    a = np.sum(y)
    #print(f"Probability of k messages (a_k):{a}")
    # a is a probability and must be < 1
    assert(a <= 1.0)

    # Probability that exactly k messages arrive during a service time
    return a


def calculate_pi(K, arrival_rate, service_distribution_pdf):
    pi_prime = np.zeros((K))
    pi = np.zeros((K))

    # Calculate the arrival (a_k) probabilities
    a = np.zeros((K))
    for k in range(0, K):
        a[k] = calculate_a_k(k, arrival_rate, service_distribution_pdf)

    #print(f"Probability of k messages arriving: {a}")
    pi_prime[0] = 1
    for k in range(0, K-1):
        # Compute every value of pi'
        s = 0
        for j in range(1, k+1):
            s += pi_prime[j] * a[k-j+1]
        pi_prime[k+1] = (1.0/a[0]) * (pi_prime[k] - s - a[k])
    #print(f"pi_prime:{pi_prime}")

    pi[0] = 1.0 / np.sum(pi_prime)
    for k in range(1, K):
        pi[k] = pi[0] * pi_prime[k]

    #print(f"pi[0]:{pi}")
    return pi


def create_pdf(samples):
    if len(samples) == 0:
        return None
        #return np.array(0), np.array(0)

    samples.sort()
    counts = []
    values = []
    prev = None
    for v in samples:
        if v == prev:
            counts[-1] += 1
        else:
            values.append(v)
            counts.append(1)
        prev = v

    values = np.array(values)
    counts = np.array(counts)
    probabilities = np.array(counts) / np.sum(counts)

    # Total probability should roughly=1
    assert abs(np.sum(probabilities) - 1.0) < 0.001
    assert values.shape == probabilities.shape

    return values, probabilities


import scipy


def calculate_pk_mmk1n(K, N, ro):
    assert(K > 0)
    assert(ro >= 0)
    p = np.zeros((K+1,))

    if (N > 0):
        print(f"\tN packets:{N} buffer size:{K}")
        for k in range(0, K+1):
            p[k] = scipy.special.binom(N, k) * (ro**k) * math.factorial(k)
        p = p / np.sum(p)
        assert(p[0] > 0)  # Normalize to 1
        p = np.clip(p, 0, 1)
        print(f":\tp0:{p[0]:e}")
        print(f"\tpk:{p}")

    return p


def calculate_queue_blocking(K, arrival_rate, mean_service_time,
                             N=None, service_pdf=None, sim_time=None):
    #print("*** Calculate queue blocking ***")
    #print(f"arrival_rate: {arrival_rate:e}")
    #print(f"mean service time: {mean_service_time:e}")
    #assert(mean_service_time > 0)
    #service_rate = 1 / mean_service_time

    ro = arrival_rate * mean_service_time
    print(f"arrival_rate:{arrival_rate:e} mean_service_time:{mean_service_time}")
    """
    if service_pdf is not None:
        # Probability density function was given, calculate parameters for an
        #  M/G/1/K queue
        # TODO: compare results with and without this
        #print(f"Calculating stats for queue of size: {K}")

        #assert(np.sum(service_pdf[1]) == 1.0)
        pi = calculate_pi(K, arrival_rate, service_pdf)

        probability_blocking = 1.0 - (1.0 / (pi[0] + ro))
        #print(f"pi0: {pi[0]} ro:{ro} probability_blocking:{probability_blocking}")

        waiting_sum = 0
        for k in range(1, K):
            waiting_sum += pi[k]
        mean_waiting_time = (waiting_sum/arrival_rate) + ((K/arrival_rate)*(pi[0] + ro - 1)) - mean_service_time
        #print(f"mean waiting time:{mean_waiting_time} mean_service_time:{mean_service_time}")
    else:
    """
    """
    if N is not None:
        print("MM1KN")
        # Calculate p0
        server_utilization = ro  # Rho in queueing theory
        # TODO: we need to reevaluate p[k] for size N-1?
        # TODO: simplify the formula, we don't need to do the factorial.
        #  I think we actually want perm nPr -- or N_Pk
        #  way
        print("Calculating Pk(N-1, 1, K)")
        p = calculate_pk_mmk1n(K, N-1, ro)
        probability_blocking = p[K]
        probability_blocking = np.nan_to_num(probability_blocking)
        #exit()
        print("Calculating Pk(N, 1, K)")
        p = calculate_pk_mmk1n(K, N, ro)
        effective_throughput = arrival_rate * (1 - probability_blocking)  # lambda_e
        #total_length = 0  # L
        queue_length = 0
        for k in range(1, K+1):
            queue_length += (k-1) * p[k]
            #total_length += k * p[k]

        total_length = queue_length
        if total_length > 0:
            total_wait = total_length / effective_throughput  # W
        else:
            total_wait = 0  # W
        mean_waiting_time = total_wait - mean_service_time  # Wq
        effective_utilization = effective_throughput * mean_service_time # rho_e
        #queue_length = total_length - effective_utilization

        print(f"mean service time:{mean_service_time} arrival rate:{arrival_rate:e} server_utilization:{server_utilization}")
        print(f"prob blocking:{probability_blocking:e} total length:{total_length} total waiting time:{total_wait} effective_throughput:{effective_throughput:e}")
        print(f"queue length:{queue_length} estimated waiting time:{mean_waiting_time}")
    """
    #else:
    print("\tMM1K")
    # No probability density function given, assume M/M/K/1 queue
    server_utilization = ro  # Rho in queueing theory
    #print(f"server utilization: {server_utilization}")
    if (server_utilization == 1):
        probability_blocking = 1 / (K+1)
        total_length = K/2
    else:
        probability_blocking = ((1 - ro)*(ro**K)) / (1 - (ro**(K+1)))
        total_length = (server_utilization/(1.0-server_utilization))
        total_length -= (((K+1) * (server_utilization**(K+1))) /
                            (1 - server_utilization**(K+1)))  # L

    #p = np.zeros((K+1))
    #p[0] = (1-ro) / (1-(ro**(K+1)))
    #for k in range(0, K+1):
    #    p[k] = ro**k * p[0]
    #print(f"p:{p}")

    if N is not None and mean_service_time != 0:
        mean_service_rate = (1/mean_service_time)
        print(f"\tservice rate:{mean_service_rate:e}")
        #time_until_saturation = K / (arrival_rate - mean_service_rate)
        #print(f"time until saturation:{time_until_saturation}")
        #time_until_all_sent = N / arrival_rate
        #print(f"N:{N} K:{K} time until all sent:{time_until_all_sent}")

        #if sim_time is not None:
        #steady_state_ratio = (time_until_all_sent - time_until_saturation) / time_until_all_sent
        #steady_state_ratio = np.clip(steady_state_ratio, 0, 1)

        # TODO: this ratio should change depending on the average depth of the link
        #  i.e., the average capacity of downstream links. The number of packets
        #  needs to be greater than this total capacity before we will see any
        #  blocking effects
        # N: Messages sent in flow
        # K: Buffer size
        steady_state_ratio2 = (N - K) / N
        steady_state_ratio2 = np.clip(steady_state_ratio2, 0, 1)
        # TODO: HACK - to disable
        #steady_state_ratio2 = 1.0
        print(f"\tsim time:{sim_time} steady state ratio2:{steady_state_ratio2}")
        #input()
        probability_blocking *= steady_state_ratio2


    # lambda_e == effective throughput
    effective_throughput = arrival_rate * (1 - probability_blocking)
    if total_length > 0:
        total_wait = total_length / effective_throughput
    else:
        total_wait = 0
    mean_waiting_time = total_wait - mean_service_time

    effective_utilization = effective_throughput * mean_service_time # rho_e
    queue_length = total_length - effective_utilization

    probability_blocking = np.clip(probability_blocking, 0.0, 1.0)
    mean_waiting_time = np.clip(mean_waiting_time, 0.0, None)
    queue_length = np.clip(queue_length, 0, None)

    print(f"\tserver_utilization:{server_utilization} total length:{total_length} total wait:{total_wait}")
    print(f"\tprob blocking {probability_blocking} queue length:{queue_length} effective_throughput:{effective_throughput:e} waiting_time:{mean_waiting_time}")
    ##print(f"mean service rate:{mean_service_rate}")

    return probability_blocking, mean_waiting_time, queue_length


def update_buffer_queue(dependencies, link, link_arrival_rates,
                        link_server_time, contention_waiting_time,
                        router_link_flows, path_arrival_rates, flows,
                        prob_link_blocking, link_arrival_count,
                        messages_buffered, max_neuron_processing):
    print("** Update buffer queue **")
    links_out = dependencies.out_edges(link)
    print(f"links_out:{links_out}")
    link_idx = reverse_graph_index(link)

    # Calculate the average server rate of all downstream links
    for downstream in links_out:
        downstream_link = downstream[1]
        print(f"downstream link:{downstream_link}")
        downstream_idx = reverse_graph_index(downstream_link)

        arrival_rate_between_links = 0
        # Find all flows through this link, going to the downstream link
        flows_in_both_links = list(set(router_link_flows[link]) &
                               set(router_link_flows[downstream_link]))
        for f in flows_in_both_links:
            path_idx = flows[f]
            arrival_rate_between_links += path_arrival_rates[path_idx[0],
                                                             path_idx[1]]
            print(f"flow:{f} path idx:{path_idx} arrival:{path_arrival_rates[path_idx[0], path_idx[1]]:e}")

        print(f"contention waiting time:{contention_waiting_time[downstream_idx]:e}")
        print(f"arrival rate between links {link_idx}->{downstream_idx} = "
              f"{arrival_rate_between_links:e}")
        print(f"total link arrival rate: {link_arrival_rates[link_idx]:e}")
        print(f"scaled contention time: {((contention_waiting_time[downstream_idx]*arrival_rate_between_links) / link_arrival_rates[link_idx]):e}")
        print(f"prob blocking:{prob_link_blocking[downstream_idx]}")
        #add = (contention_waiting_time[downstream_idx] + \
        #    (1 / (1 - prob_link_blocking[downstream_idx]))) * \
        #    arrival_rate_between_links
        #link_server_time[link_idx] += \
        #    (contention_waiting_time[downstream_idx] + \
        #    (1 / (1 - prob_link_blocking[downstream_idx]))) * \
        #    arrival_rate_between_links
        link_server_time[link_idx] += \
            (arrival_rate_between_links) * contention_waiting_time[downstream_idx]

    # Normalize server time with respect with the total flow going through this
    # link
    if (len(links_out) > 0) and (link_arrival_rates[link_idx] > 0):
        link_server_time[link_idx] /= link_arrival_rates[link_idx]
    print(f"buffer server time: {link_server_time[link_idx]}")
    assert(link_server_time[link_idx] != math.nan)
    assert(link_server_time[link_idx] < 1.0)
    # else this is a leaf link, and the service time is already given by the path

    #"""
    if link_idx[2] == 0 or link_idx[2] == 2:  # Network tile to network tile
        buffer_size = 16
    elif link_idx[2] == 1 or link_idx[2] == 3:
        buffer_size = 10
    elif (link_idx[2] >= 4 and link_idx[2] < 8):  # Cores to network
        buffer_size = 8
    elif link_idx[2] >= 8:  # Network to cores
        buffer_size = 24
    #"""

    prob_blocked, link_waiting_time, queue_length = \
        calculate_queue_blocking(buffer_size, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None,
                                 sim_time=max_neuron_processing)
    #prob_blocked, link_waiting_time, queue_length = \
    #    calculate_queue_blocking(buffer_size, link_arrival_rates[link_idx],
    #                             link_server_time[link_idx],
    #                             N=None,
    #                             service_pdf=None,
    #                             sim_time=max_neuron_processing)

    """
    prob_blocked1, link_waiting_time1 = \
        calculate_queue_blocking(10, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None)
    prob_blocked2, link_waiting_time2 = \
        calculate_queue_blocking(100, link_arrival_rates[link_idx],
                                 link_server_time[link_idx],
                                 N=link_arrival_count[link_idx],
                                 service_pdf=None)
    print(f"prob block 1:{prob_blocked1} 2:{prob_blocked2}")
    """
    #input()

    print(f"link buffer waiting:{link_waiting_time:e}")
    prob_link_blocking[link_idx] = np.clip(prob_blocked, 0, 1)
    messages_buffered[link_idx] = np.clip(queue_length, 0, None)
    return link_waiting_time


def update_contention_queue(dependencies, link, link_arrival_rates,
                            contention_waiting_time,
                            prob_link_blocking,
                            link_server_time,
                            sim_time):
    print("** Update contention queue **")
    links_in = dependencies.in_edges(link)
    print(f"links_in:{links_in}")
    link_in_count = len(links_in)
    print(f"link in count:{link_in_count}")
    link_idx = reverse_graph_index(link)

    # The server time of the queue is just the average time that the link
    #  is blocked for
    time_blocked = 0
    if (link_in_count >= 1) and link_arrival_rates[link_idx] > 0:
        # This equation doesn't work / make sense
        #time_blocked = (1 / link_arrival_rates[link_idx]) * (1 / (1 - prob_link_blocking[link_idx]))
        time_blocked = prob_link_blocking[link_idx] * link_server_time[link_idx]
        print(f"time_blocked:{time_blocked} arrival:{link_arrival_rates[link_idx]:e} prob_block:{prob_link_blocking[link_idx]}")
        assert(time_blocked != math.nan)
        assert(time_blocked >= 0)
        contention_server_time = time_blocked
        print(f"prob blocked: {prob_link_blocking[link_idx]} time blocked:{contention_server_time:e}")

        _, contention_waiting_time[link_idx], _ = \
            calculate_queue_blocking(link_in_count, link_arrival_rates[link_idx],
                                    contention_server_time, None,
                                    sim_time=sim_time)
    print(f"contention waiting time:{contention_waiting_time[link_idx]}")
    print ("**end of contention**")

    return


n_tiles = 32
router_link_names = ("north", "east", "south", "west",
                    "core_1_to_net", "core_2_to_net",
                    "core_3_to_net", "core_4_to_net",
                    "net_to_core_1", "net_to_core_2",
                    "net_to_core_3", "net_to_core_4")


def sim_delay_hops(flows):
    """Calculate delay of each message based on the hop count

    The simplest delay calculation scheme, beyond just setting network delay=0
    """
    delays = np.zeros((128, 128))
    cost_per_hop = 4.1e-9  # s
    for src_core, dest_core in flows:
        src_tile = src_core // 4
        dest_tile = dest_core // 4
        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4
        x_hops = abs(src_x - dest_x)
        y_hops = abs(src_y - dest_y)
        delays[src_core][dest_core] = (x_hops + y_hops) * cost_per_hop

    return delays

import enum
def sim_schedule_event_based(messages):
    # Model the queues in each router link as a set of FIFOs. Here we will
    #  essentially schedule every single event as messages go from router to
    #  router.
    noc_buffers = []
    for x in range(0, 8):
        y_buffers = []
        for y in range(0, 4):
            link_buffers = [[] for _ in range(0, 12)]
            y_buffers.append(link_buffers)
        noc_buffers.append(y_buffers)

    class State(enum.Enum):
        SEND = 1
        RECEIVE = 2

    def get_next_link(m):
        if m.pos is None:
            # if pos is None, the next link is the first one
            src_tile = m.src_core // 4
            x = src_tile // 4
            y = src_tile % 4
            link = 4 + m.src_core % 4
        else:
            x, y, link = m.pos
            dest_tile = m.dest_core // 4
            dest_x = dest_tile // 4
            dest_y = dest_tile % 4
            if x < dest_x:
                x = x + 1
            elif x > dest_x:
                x = x - 1
            elif y < dest_y:
                y = y + 1
            elif y > dest_y:
                y = y - 1

            # Work out the direction of the next link
            if x < dest_x:
                link = 1
            elif x > dest_x:
                link = 3
            elif y < dest_y:
                link = 0
            elif y > dest_y:
                link = 2
            else:
                link = 8 + m.dest_core % 4

        # If the next link is full, it will block us from sending
        buffer_sizes = (16, 10, 16, 10, 8, 8, 8, 8, 24, 24, 24, 24)
        if len(noc_buffers[x][y][link]) >= buffer_sizes[link]:
            return None  # to indicate the next link is full
        else:
            return x, y, link


    def receive_message(m, t):
        # Schedule an event to update the noc once we finish processing
        print(f"Creating a message event for {m} at time {t}")
        E = {"message": m, "type": State.RECEIVE}
        priority.append((E, t + m.receive_delay))
        priority.sort(key=lambda x: x[1], reverse=True)
        return


    def send_message(m, t):
        print(f"Sending message {m} from cid:{m.src_core}->cid:{m.dest_core}")
        dest_pos = get_next_link(m)
        if dest_pos is None:
            # Message is being blocked by NoC, we need to push this
            return False

        while (dest_pos is not None and
               len(noc_buffers[dest_pos[0]][dest_pos[1]][dest_pos[2]]) == 0 and
               dest_pos[2] < 8):
            m.pos = dest_pos
            dest_pos = get_next_link(m)

        dest_x, dest_y, dest_link = dest_pos
        if dest_link >= 8 and len(noc_buffers[dest_x][dest_y][dest_link]) == 0:
            # We successfully sent the message across the NoC, no blocking
            #  We can start receiving the message
            receive_message(m, t)
        else:
            print(f"Message put in queue noc[{dest_x}][{dest_y}][{dest_link}] "
                  f"Queue len:{len(noc_buffers[dest_x][dest_y][dest_link])}")
            noc_buffers[dest_x][dest_y][dest_link].append(m)

        print(f"Message sent to {dest_pos}")

    def update_message(m):
        if src_pos is None:
            pass
        else:
            dest_x, dest_y, dest_link = src_pos

        next_x, next_y, next_link = next_pos
        dest_pos = src_pos
        src_x, src_y, src_link = src_pos
        dest_x, dest_y, dest_link = dest_pos
        while (len(noc_buffers[dest_x][dest_y][dest_link]) == 0 and
               next_pos is not None):
            dest_x, dest_y, dest_link = next_x, next_y, next_link
            m.pos = next_x, next_y, next_link
            next_x, next_y, next_link = get_next_link(m)

        print(f"m.pos: {m.pos}->{(dest_x, dest_y, dest_link)}")

        if dest_link > 8 and noc_buffers[dest_x][dest_y][dest_link] == 0:
            # We successfully sent the message across the NoC, no blocking
            #  We can start receiving the message
            receive_message(m)
        elif m.pos != src_pos:
            # The message made it some way through the NoC but not to the end,
            #  so update the corresponding queues
            if src_pos is not None:
                noc_buffers[src_x][src_y][src_link].pop()
            noc_buffers[dest_x][dest_y][dest_link].append()

        exit()
        return True

    def schedule_next_message(m, t):
        # Schedule the next message
        next_message = messages[m.src_core][-1]
        print(f"Scheduling next message for core:{m.src_core}")
        E = {"message": next_message, "type": State.SEND}
        priority.append((E, t + next_message.generation_delay))
        priority.sort(key=lambda x: x[1], reverse=True)
        return

    # Initialize priority queue by pushing events to send the first messages
    priority = []
    for core in range(0, 128):
        if len(messages[core]) > 0:
            first_message = messages[core][-1]
            E = {"message": first_message, "type": State.SEND}
            priority.append((E, first_message.generation_delay))
    priority.sort(key=lambda x: x[1], reverse=True)

    while len(priority) > 0:
        E, t = priority.pop()
        print(E)
        m = E["message"]

        if E["type"] == State.SEND:
            # First, figure out if we are able to send a message or not. If the NoC
            #  is blocking, put the message into a pending state
            m = messages[m.src_core].pop()
            assert(m == E["message"])
            print(f"src core:{m.src_core} pos:{m.pos}")
            m.pending = send_message(m, t)
            # Message has been sent successfully, schedule the next message
            if not m.pending and len(messages[m.src_core]) > 0:
                schedule_next_message(m, t)
        elif E["type"] == State.RECEIVE:  # State.Receive
            # Update the state of the NoC. While we encounter queues, pop an
            #  element from them. If the queue is full, also go to its
            #  connected neighbours and pop.
            pass
    print("*** Finished modeling NoC ***")
    exit()
    return


def sim_delay_mm1k(df):
    n_cores = 128
    path_counts = np.zeros((n_cores, n_cores), dtype=int)  # [src core, dest core]
    path_arrival_latencies = np.zeros((n_cores, n_cores))
    path_server_mean_latencies = np.zeros((n_cores, n_cores))
    router_link_flows = [[] for _ in range(0, 8*4*12)]

    path_server_latencies = [[] for _ in range(0, 128)]
    path_server_pdf = [[] for _ in range(0, 128)]
    for s in range(0, 128):
        path_server_latencies[s] = [[] for _ in range(0, 128)]
        path_server_pdf[s] = [[] for _ in range(0, 128)]

    router_link_counts = np.zeros((8, 4, 12))
    router_link_arrival_rates = np.zeros((8, 4, 12))

    dependencies = nx.DiGraph()
    # Flatten graph index into 1-dimension
    dependencies.add_nodes_from(range(0, 8*4*12))

    # Parse the rate of messages between two cores in the network and build a
    #  dependency graph
    for _, row in df.iterrows():
        src_tile, src_core = hw_str_to_core(row["src_hw"])
        dest_tile, dest_core = hw_str_to_core(row["dest_hw"])
        path_arrival_latencies[src_core, dest_core] += row["generation_delay"]
        #if row["processing_latency" != 0]:
        #    print(f"latency:{path_arrival_latencies[src_core, dest_core]}")
        assert(row["processing_latency"] >= 0)
        assert(row["processing_latency"] != math.nan)
        path_server_mean_latencies[src_core, dest_core] += \
            row["processing_latency"]
        assert(np.all(path_server_mean_latencies >= 0))
        path_counts[src_core, dest_core] += 1
        path_server_latencies[src_core][dest_core].append(row["processing_latency"])

    #assert(np.all(path_server_mean_latencies >= 0))
    #assert(not np.any(np.isnan(path_server_mean_latencies)))
    #assert(not np.any(np.isnan(path_counts)))
    #assert(np.all(path_counts >= 0))

    path_server_mean_latencies = np.divide(path_server_mean_latencies,
                                           path_counts, where=(path_counts>=1),
                                           dtype=float)
    np.clip(path_server_mean_latencies, 0.0, None)
    path_server_mean_latencies = np.nan_to_num(path_server_mean_latencies)

    #path_arrival_rates = np.divide(path_counts, np.max(path_arrival_latencies),
    #                               where=path_counts>0)  # packets/s
    path_arrival_rates = np.divide(path_counts, path_arrival_latencies*10,
                                   where=path_counts>0)  # packets/s

    # Treat the neuron processing time as the simulation window
    max_neuron_processing = np.max(np.sum(path_arrival_latencies, axis=1))

    #print(np.ma.masked_less(path_server_rates, 1.0))
    flows = np.argwhere(path_counts >= 1)
    # Figure out the distribution of latencies for each flow
    for s in range(0, 128):
        for d in range(0, 128):
            path_server_pdf[s][d] = create_pdf(path_server_latencies[s][d])

    #print(path_rates)

    # *** Router link model variables ***
    prob_link_blocking = np.zeros((8, 4, 12))
    mean_link_waiting_time = np.zeros((8, 4, 12))
    mean_link_service_time = np.zeros((8, 4, 12))
    messages_buffered = np.zeros((8, 4, 12))

    # *** Contention model variables ***
    contention_waiting_time = np.zeros((8, 4, 12))

    # Now parse all the links between two cores. Every router has 12 links in
    #  the NoC: 4 links in the going North, East, south and West. Every router
    #  has 4 links going from the core to the router, and 4 links going out to
    #  the cores. Compute the packet arrival rates for all links on all paths
    for idx, flow in enumerate(flows):
        src_core, dest_core = flow
        assert(src_core < 128)
        assert(dest_core < 128)

        src_tile = src_core // 4
        dest_tile = dest_core // 4
        assert(src_tile < 32)
        assert(dest_tile < 32)

        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4
        path_rate = path_arrival_rates[src_core, dest_core]
        path_count = path_counts[src_core, dest_core]

        # Account for all the hops in between the src and dest tile
        initial_direction = 4 + (src_core % 4)
        router_link_counts[src_x, src_y, initial_direction] += path_count
        router_link_arrival_rates[src_x, src_y, initial_direction] += path_rate
        router_link_flows[graph_index(src_x, src_y, initial_direction)].append(idx)
        prev_link = (src_x, src_y, initial_direction)

        while src_x < dest_x:
            src_x += 1
            router_link_counts[src_x, src_y, 1] += path_count  # east
            router_link_arrival_rates[src_x, src_y, 1] += path_rate
            router_link_flows[graph_index(src_x, src_y, 1)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 1))
            prev_link = (src_x, src_y, 1)

        while src_x > dest_x:
            src_x -= 1
            router_link_counts[src_x, src_y, 3] += path_count  # west
            router_link_arrival_rates[src_x, src_y, 3] += path_rate
            router_link_flows[graph_index(src_x, src_y, 3)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 3))
            prev_link = (src_x, src_y, 3)

        while src_y < dest_y:
            src_y += 1
            router_link_counts[src_x, src_y, 0] += path_count  # north
            router_link_arrival_rates[src_x, src_y, 0] += path_rate
            router_link_flows[graph_index(src_x, src_y, 0)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 0))
            prev_link = (src_x, src_y, 0)

        while src_y > dest_y:
            src_y -= 1
            router_link_counts[src_x, src_y, 2] += path_count  # south
            router_link_arrival_rates[src_x, src_y, 2] += path_rate
            router_link_flows[graph_index(src_x, src_y, 2)].append(idx)
            dependencies.add_edge(graph_index(*prev_link),
                                graph_index(src_x, src_y, 2))
            prev_link = (src_x, src_y, 2)

        final_direction = 8 + (dest_core % 4)
        router_link_counts[src_x, src_y, final_direction] += path_count
        router_link_arrival_rates[src_x, src_y, final_direction] += path_rate
        router_link_flows[graph_index(src_x, src_y, final_direction)].append(idx)
        dependencies.add_edge(graph_index(*prev_link),
                            graph_index(src_x, src_y, final_direction))

        # Track the service time at the receiving link
        mean_link_service_time[src_x, src_y, final_direction] = \
            path_server_mean_latencies[src_core, dest_core]
        #print(mean_contention_delay[src_x, src_y, final_direction])


    # Create dependency graph of all links
    sorted_links = list(nx.topological_sort(dependencies))
    sorted_link_info = [reverse_graph_index(n) for n in sorted_links]
    nx.nx_agraph.write_dot(dependencies, "runs/dependencies.dot")
    D = nx.nx_agraph.to_agraph(dependencies)
    # If we want to plot dependencies
    #plt.figure()
    #pos=graphviz_layout(dependencies, prog='dot')
    #nx.draw(dependencies, pos, with_labels=False, arrows=True, node_size=0.5)


    # Now iterate over all flows again in reverse topological order, and
    #  calculate the contention delay followed by the service time at each
    #  link queue. Calculate the total arrival rate at each queue.
        # Calculate the average contention delay for this link


    sorted_links.reverse()
    #print(sorted_links)
    for link in sorted_links:
        print(f"UPDATING LINK: {link}")
        link_wait_time = \
            update_buffer_queue(dependencies, link, router_link_arrival_rates,
                                mean_link_service_time,
                                contention_waiting_time,
                                router_link_flows,
                                path_arrival_rates, flows,
                                prob_link_blocking,
                                router_link_counts,
                                messages_buffered,
                                max_neuron_processing)
        mean_link_waiting_time[reverse_graph_index(link)] = link_wait_time
        update_contention_queue(dependencies, link, router_link_arrival_rates,
                                contention_waiting_time, prob_link_blocking,
                                mean_link_service_time,
                                max_neuron_processing)


    # Now go through all flows, and calculate the delay for each flows path
    flow_latencies = np.zeros((128, 128))
    for flow in flows:
        src_core, dest_core = flow
        src_tile = src_core // 4
        dest_tile = dest_core // 4
        assert(src_tile < 32)
        assert(dest_tile < 32)

        src_x = src_tile // 4
        src_y = src_tile % 4
        dest_x = dest_tile // 4
        dest_y = dest_tile % 4

            # Account for all the hops in between the src and dest tile
        initial_direction = 4 + (src_core % 4)
        flow_latencies[src_core, dest_core] += \
            mean_link_waiting_time[src_x, src_y, initial_direction]

        while src_x < dest_x:
            src_x += 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 1]

        while src_x > dest_x:
            src_x -= 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 3]

        while src_y < dest_y:
            src_y += 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 0]

        while src_y > dest_y:
            src_y -= 1
            flow_latencies[src_core, dest_core] += \
                mean_link_waiting_time[src_x, src_y, 2]

        final_direction = 8 + (dest_core % 4)
        flow_latencies[src_core, dest_core] += \
            mean_link_waiting_time[dest_x, dest_y, final_direction]

    """
    # Plot a heat map
    plt.figure()
    plt.title("Path Arrival Rates")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=path_arrival_rates, cmap="YlOrRd", s=0.4)
    plt.colorbar()

    import matplotlib.ticker as ticker
    def create_subplots(c, title=""):
        fig, ax = plt.subplots(nrows=3, ncols=4)
        cmin = 0
        cmax = np.max(c)
        plt.suptitle(title)

        for i in range(0, 12):
            x, y = np.meshgrid(np.linspace(0, 7, 8, dtype=int),
                            np.linspace(0, 3, 4, dtype=int))
            pcm = ax[i//4,i%4].scatter(x, y, c=c[:,:,i], cmap="YlOrRd",
                                       vmin=cmin, vmax=cmax)
            ax[i//4,i%4].set_title(router_link_names[i])
            ax[i//4,i%4].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

            fig.colorbar(pcm, ax=ax[i//4,i%4])

    #print(router_link_counts[:,:,4])
    #print(router_link_arrival_latency)
    #print(router_link_rates)
    create_subplots(router_link_arrival_rates, "Link Arrival Rates")
    create_subplots(prob_link_blocking, "Probability of Blocking")
    create_subplots(mean_link_waiting_time, "Mean Link Wait Time")
    create_subplots(messages_buffered, "Messaged Buffered")
    create_subplots(router_link_counts, "Link Counts")

    # Plot a heat map
    plt.figure()
    plt.title("Flow Latencies")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=flow_latencies, cmap="YlOrRd", s=0.4)
    plt.colorbar()

    plt.figure()
    plt.title("Path Mean Receive Latencies")
    plt.xlabel("Source Core")
    plt.ylabel("Destination Core")
    x, y = np.meshgrid(np.linspace(0, 127, 128), np.linspace(0, 127, 128))
    plt.scatter(x, y, c=path_server_mean_latencies, cmap="YlOrRd", s=0.4)
    plt.colorbar()
    """

    #print(f"mean link transfer delay: {flow_latencies}")
    return flow_latencies, path_counts, path_server_mean_latencies


def schedule_messages_simple(messages):
    neuron_processing = np.zeros((128,))
    message_receiving = np.zeros((128,))
    t = 0.0
    for m in messages:
        neuron_processing[m.src_core] += m.generation_delay
        message_receiving[m.dest_core] += m.receive_delay

    for i in range(0, 128):
        t = max((neuron_processing[i], message_receiving[i], t))
    return t


def schedule_messages_detailed(messages):
    priority = []
    receiving_pipeline_busy = np.zeros((128,))
    # Setup the priority queue, using just a normal list that will be sorted
    #  after every insert. This was easier than figuring out the PriorityQueue
    #  class for floating point priority...
    t = 0.0
    for core in range(0, 128):
        if len(messages[core]) > 0:
            first_message = messages[core][-1]
            priority.append((core, first_message.generation_delay))
    priority.sort(key=lambda x: x[1], reverse=True)

    while len(priority) > 0:
        core, t = priority.pop()
        m = messages[core].pop()
        earliest_receive_time = max(t + m.network_delay,
                                    receiving_pipeline_busy[m.dest_core])
        receiving_pipeline_busy[m.dest_core] = \
            earliest_receive_time + m.receive_delay

        # Push next message
        if len(messages[core]) > 0:
            next_message = messages[core][-1]
            priority.append((core, t+next_message.generation_delay))
            priority.sort(key=lambda x: x[1], reverse=True)

    max_receive = np.max(receiving_pipeline_busy)
    return max(t, max_receive)


class Message:
    """Spike message timing information"""
    def __init__(self, generation_delay=None, receive_delay=None,
                 src_hw=None, dest_hw=None, hops=None):
        self.generation_delay = float(generation_delay)
        self.receive_delay = float(receive_delay)
        self.src_tile, self.src_core = hw_str_to_core(src_hw)
        self.dest_tile, self.dest_core = hw_str_to_core(dest_hw)
        self.network_delay = None
        self.hops = hops
        self.pos = None
        return

    def __repr__(self):
        return str(self) + "\n"

    def __str__(self):
        s = ""
        for k,v in self.__dict__.items():
            s += f"{k}:{v} "
        return s


# 1. Read in the network
filename = "dvs_messages.trace"
#filename = "latin_messages.trace"
df = pd.read_csv(filename,
                 converters={"src_hw": str, "dest_hw": str})

timesteps = 128
max_latencies = np.zeros((timesteps,))
mean_latencies = np.zeros((timesteps,))
total_flow_latencies = np.zeros((timesteps,))
total_core_latencies = np.zeros((timesteps,))
max_synapse_processing = np.zeros((timesteps,))
max_neuron_processing = np.zeros((timesteps,))
scheduled_latency = np.zeros((timesteps,))
flow_delays1 = np.zeros((timesteps,))
path_counts = np.zeros((128, 128), dtype=int)  # [src core, dest core]

# Queue of messages for each core
messages = [[] for _ in range(0, 128)]

#for timestep in range(0, timesteps):
for timestep in range(102, 103):
    message_generation_latencies = np.zeros((128, 128))
    message_receive_latencies = np.zeros((128, 128))
    df_timestep = df[df["timestep"] == timestep]

    for i, row in df_timestep.iterrows():
        src_tile, src_core = hw_str_to_core(row["src_hw"])
        dest_tile, dest_core = hw_str_to_core(row["dest_hw"])
        messages[src_core].append(Message(row["generation_delay"],
                                       row["processing_latency"],
                                       row["src_hw"], row["dest_hw"],
                                       row["hops"]))

        #message_generation_latencies[src_core, dest_core] += row["generation_delay"]
        if message_generation_latencies[src_core, dest_core] == 0:
            message_generation_latencies[src_core, dest_core] = \
                row["generation_delay"]
        else:
            message_generation_latencies[src_core, dest_core] += \
                min(message_generation_latencies[src_core, dest_core],
                    row["generation_delay"])
        message_receive_latencies[src_core, dest_core] += row["processing_latency"]
        path_counts[src_core, dest_core] += 1

    print(f"Scheduling messages for timestep:{timestep}")
    sim_schedule_event_based(messages)

    #TODO explore
    #message_generation_latencies = df_timestep["generation_delay"].min()
    #message_generation_latencies = df_timestep["generation_delay"].min()

    flows = np.argwhere(path_counts >= 1)
    flow_delays1[timestep] = np.max(sim_delay_hops(flows))
    flow_delays, flow_counts, _ = sim_delay_mm1k(df_timestep)

    # For all messages, update their network delay based on the average flow
    #  delay
    for core in range(0, 128):
        for m in messages[core]:
            m.network_delay = flow_delays[core, m.dest_core]
    scheduled_latency[timestep] = schedule_messages_detailed(messages)

    print(f"i:{timestep} max flow delay: {np.max(flow_delays):e}")
    max_latencies[timestep] = np.max(flow_delays)

    #print(f"mean server:{message_receive_latencies}")
    #print(f"counts:{flow_counts}")

    max_neuron_processing[timestep] = np.max(np.sum(message_generation_latencies, axis=1))
    max_synapse_processing[timestep] = np.max(np.sum(message_receive_latencies, axis=0))


PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
LOIHI_TIME_DATA_FILENAME = "loihi_gesture_32x32_time.csv"
NETWORK_DIR = os.path.join(PROJECT_DIR, "runs", "dvs", "loihi_gesture_32x32")
DVS_RUN_DIR = os.path.join(PROJECT_DIR, "runs", "dvs")
LOIHI_TIME_DATA_PATH = os.path.join(DVS_RUN_DIR, LOIHI_TIME_DATA_FILENAME)

plt.figure()
loihi_data = pd.read_csv(LOIHI_TIME_DATA_PATH)
loihi_times = np.array(loihi_data.loc[:, :] / 1.0e6)
plt.plot(np.arange(1, timesteps-1), loihi_times[0:(timesteps-2), 0] * 1.0e6, "-")
plt.plot(max_latencies[1:] * 1.0e6)
plt.plot(max_synapse_processing[1:] * 1.0e6, "--")
plt.plot(max_neuron_processing[1:] * 1.0e6, "--")
plt.plot(scheduled_latency[1:] * 1.0e6)
plt.plot(flow_delays1[1:] * 1.0e6)

#plt.legend(("Measured", "Max", "Mean", "Max Synapse", "Max Neuron"), fontsize=7)
plt.legend(("Measured", "Max Path Delay", "Max. Message Receiving", "Max. Neuron Processing", "Scheduled", "Simple hop latency"), fontsize=7)
plt.ylabel("Time-step Latency ($\mu$s)")
plt.xlabel("Time-step")
plt.yticks(np.arange(0, 61, 10))

plt.show()
