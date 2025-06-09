// Copyright (c) 2025 - The University of Texas at Austin
//  This work was produced under contract #2317831 to National Technology and
//  Engineering Solutions of Sandia, LLC which is under contract
//  No. DE-NA0003525 with the U.S. Department of Energy.
// schedule.hpp: Schedule global order of messages on a neuromorphic chip

#include <memory>

#include "message.hpp"
#include "timestep.hpp"

sanafe::Timestep::Timestep(const long int ts)
        : timestep(ts)
{
}

void sanafe::Timestep::set_cores(const size_t core_count)
{
    messages.resize(core_count);
}

sanafe::TimestepHandle::TimestepHandle()
{
    // Initialize this empty object without any data. The assumption is that an
    //  empty object will be assigned to something else before its used.
}

sanafe::TimestepHandle::TimestepHandle(const long int timestep_num)
        : handle(std::make_shared<Timestep>(timestep_num))
{
}