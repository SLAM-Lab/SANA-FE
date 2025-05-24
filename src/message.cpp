#include <cstddef>

#include "chip.hpp"
#include "core.hpp"
#include "mapped.hpp"
#include "message.hpp"
#include "tile.hpp"

sanafe::Message::Message(const long int id, const SpikingChip &hw,
        const MappedNeuron &n, const long int timestep)
        : timestep(timestep)
        , mid(id)
        , src_neuron_id(n.id)
        , src_neuron_group_id(n.parent_group_name)
{
    // If no axon was given create a message with no destination. By
    //  default, messages without destinations act as a placeholder for neuron
    //  processing
    const Core &src_core = *(n.core);
    const Tile &src_tile = hw.tiles[src_core.parent_tile_id];
    src_x = src_tile.x;
    src_y = src_tile.y;
    src_tile_id = src_tile.id;
    src_core_id = src_core.id;
    src_core_offset = src_core.offset;
}

sanafe::Message::Message(const long int id, const SpikingChip &hw,
        const size_t axon_address, const MappedNeuron &n,
        const long int timestep)
        : Message(id, hw, n, timestep)
{
    const Core &src_core = *(n.core);
    const AxonOutModel &src_axon = src_core.axons_out[axon_address];
    const Tile &dest_tile = hw.tiles[src_axon.dest_tile_id];
    const Core &dest_core = dest_tile.cores[src_axon.dest_core_offset];
    const AxonInModel &dest_axon = dest_core.axons_in[src_axon.dest_axon_id];

    placeholder = false;
    spikes = dest_axon.synapse_addresses.size();
    dest_x = dest_tile.x;
    dest_y = dest_tile.y;
    dest_tile_id = dest_tile.id;
    dest_core_id = dest_core.id;
    dest_core_offset = dest_core.offset;
    dest_axon_id = src_axon.dest_axon_id;
    dest_axon_hw = 0;
}
