// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// pytrace.hpp
#ifndef PYTRACE_HEADER_INCLUDED_
#define PYTRACE_HEADER_INCLUDED_

#include <fstream>
#include <map>
#include <variant>

#include <pybind11/pybind11.h>

#include "fwd.hpp"
#include "message.hpp"
#include "network.hpp"
#include "timestep.hpp"
#include "schedule.hpp"

// Forward declarations
using PerfStatistic = std::variant<size_t, long int, double>;

pybind11::dict message_to_dict(const sanafe::Message &m);
std::map<std::string, PerfStatistic> timestep_data_to_map(
        const sanafe::Timestep &ts);

// Determine spike trace mode
enum TraceMode : uint8_t
{
    trace_none = 0UL, // No tracing
    trace_file = 1UL, // Write to file
    trace_memory = 2UL, // Store in memory for return
};

class PyTrace
{
public:
    PyTrace(sanafe::SpikingChip *chip, pybind11::object &trace_obj,
            const bool overwrite_trace)
            : parent_chip(chip)
    {
        open(trace_obj, overwrite_trace);
    }
    virtual ~PyTrace()
    {
        if (file.is_open())
        {
            file.close();
        }
    }
    PyTrace(const PyTrace &&move) = delete;
    PyTrace &operator=(const PyTrace &copy) = delete;
    PyTrace &operator=(PyTrace &&move) = delete;

    void write_header();
    void open(pybind11::object trace_obj, const bool overwrite_trace);

    virtual pybind11::object get_python_object() const = 0;

protected:
    sanafe::SpikingChip *parent_chip;
    pybind11::object obj{pybind11::none()};
    TraceMode mode{TraceMode::trace_none};
    std::ofstream file;

private:
    virtual void get_header_string(std::ostringstream &ss) = 0;
};

class PyHardwareTrace : public PyTrace
{
public:
    PyHardwareTrace(sanafe::SpikingChip *chip, pybind11::object &trace_obj,
            const bool overwrite_trace) : PyTrace(chip, trace_obj, overwrite_trace) {}
    void record_hw_metrics(const sanafe::Timestep &ts);

private:
    virtual void write_hw_trace_to_stream(
            std::ostringstream &ss, const sanafe::Timestep &ts) = 0;
    virtual void append_hw_metrics_to_memory(const sanafe::Timestep &ts) = 0;
};

class PyNetworkTrace : public PyTrace
{
public:
    PyNetworkTrace(sanafe::SpikingChip *chip, pybind11::object &trace_obj,
            const bool overwrite_trace)
            : PyTrace(chip, trace_obj, overwrite_trace) {}
    void record_net_activity(long int timestep);

private:
    virtual void write_net_trace_to_stream(std::ostringstream &ss, long int timestep) = 0;
    virtual void append_net_activity_to_memory() = 0;
};

class PySpikeTrace : public PyNetworkTrace
{
public:
    PySpikeTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyNetworkTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PySpikeTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        sanafe::SpikingChip::sim_trace_write_spike_header(ss);
    }
    void write_net_trace_to_stream(
            std::ostringstream &ss, const long int timestep) override
    {
        parent_chip->sim_trace_record_spikes(ss, timestep);
    }
    void append_net_activity_to_memory() override
    {
        auto spikes = parent_chip->get_spikes();
        data.emplace_back(std::move(spikes));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            pybind11::list data_list;
            for (const auto &spikes : data)
            {
                pybind11::list spike_list;
                for (const sanafe::NeuronAddress &spike : spikes)
                {
                    spike_list.append(pybind11::cast(spike));
                }
                data_list.append(spike_list);
            }
            return data_list;
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<sanafe::NeuronAddress>> data;
};

class PyPotentialTrace : public PyNetworkTrace
{
public:
    PyPotentialTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyNetworkTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PyPotentialTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        parent_chip->sim_trace_write_potential_header(ss);
    }
    void write_net_trace_to_stream(
            std::ostringstream &ss, const long int timestep) override
    {
        parent_chip->sim_trace_record_potentials(ss, timestep);
    }
    void append_net_activity_to_memory() override
    {
        auto potentials = parent_chip->get_potentials();
        data.emplace_back(std::move(potentials));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            return pybind11::cast(data);
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<double>> data;
};

class PyPerfTrace : public PyHardwareTrace
{
public:
    ~PyPerfTrace() override = default;
    PyPerfTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyHardwareTrace(chip, trace_obj, overwrite_trace)
    {
    }
    void get_header_string(std::ostringstream &ss) override
    {
        parent_chip->sim_trace_write_perf_header(ss);
    }
    void write_hw_trace_to_stream(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        parent_chip->sim_trace_record_perf(ss, ts);
    }
    void append_hw_metrics_to_memory(const sanafe::Timestep &ts) override
    {
        std::map<std::string, PerfStatistic> stats = timestep_data_to_map(ts);

        // Get any optional traces and cast double values to Python objects
        const std::map<std::string, double> optional_perf_traces =
                parent_chip->sim_trace_get_optional_traces();
        for (const auto &[name, value] : optional_perf_traces)
        {
            stats[name] = value;
        }

        for (const auto &[key, value] : stats)
        {
            data[key].emplace_back(value);
        }
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            return pybind11::cast(data);
        }
        return pybind11::none();
    }

private:
    std::map<std::string, std::vector<PerfStatistic>> data;
};

class PyMessageTrace : public PyHardwareTrace
{
public:
    PyMessageTrace(sanafe::SpikingChip *chip, pybind11::object trace_obj,
            const bool overwrite_trace)
            : PyHardwareTrace(chip, trace_obj, overwrite_trace)
    {
    }
    ~PyMessageTrace() override = default;
    void get_header_string(std::ostringstream &ss) override
    {
        sanafe::SpikingChip::sim_trace_write_message_header(ss);
    }
    void write_hw_trace_to_stream(
            std::ostringstream &ss, const sanafe::Timestep &ts) override
    {
        std::vector<std::reference_wrapper<const sanafe::Message>> all_messages;
        for (const sanafe::MessageFifo &q : ts.messages)
        {
            for (const sanafe::Message &m : q)
            {
                all_messages.emplace_back(m);
            }
        }
        // Sort messages in message ID order (with placeholders last)
        std::sort(all_messages.begin(), all_messages.end(),
                sanafe::CompareMessagesByID{});
        // Save the messages in sorted order
        for (const sanafe::Message &m : all_messages)
        {
            sim_trace_record_message(ss, m);
        }
    }
    void append_hw_metrics_to_memory(const sanafe::Timestep &ts) override
    {
        std::vector<sanafe::Message> timestep_messages;

        // Not crucial, but its nice to print messages in ID order.
        //  Copy all the messages into a single vector and sort
        for (const sanafe::MessageFifo &q : ts.messages)
        {
            for (const sanafe::Message &m : q)
            {
                timestep_messages.emplace_back(m);
            }
        }
        // Sort messages in message order, starting at lowest mid first
        std::sort(timestep_messages.begin(), timestep_messages.end(),
                [](const sanafe::Message &left, const sanafe::Message &right) {
                    return left.mid < right.mid;
                });
        data.emplace_back(std::move(timestep_messages));
    }
    pybind11::object get_python_object() const override
    {
        // Acquire GIL to be safe, even though simulations should have finished
        //  at this point and the GIL should be reacquired already
        const pybind11::gil_scoped_acquire acquire;
        if (mode == TraceMode::trace_memory)
        {
            // Convert each message to a dict, meaning we need to manually
            //  iterate over every message in data
            pybind11::list data_list;
            for (const auto &ts_messages : data)
            {
                pybind11::list ts_list;
                for (const auto &m : ts_messages)
                {
                    ts_list.append(message_to_dict(m));
                }
                data_list.append(ts_list);
            }
            return data_list;
        }
        return pybind11::none();
    }

private:
    std::vector<std::vector<sanafe::Message>> data;
};

#endif