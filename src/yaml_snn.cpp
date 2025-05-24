#include <ryml.hpp> // NOLINT(misc-include-cleaner)
#include <ryml_std.hpp> // NOLINT(misc-include-cleaner)

#include "description.hpp"
#include "yaml_snn.hpp"

sanafe::SpikingNetwork sanafe::description_parse_network_file_yaml(
        std::ifstream &fp, Architecture &arch)
{
    if (!fp.is_open())
    {
        throw std::runtime_error("Error opening file\n");
    }

    // Get file size
    fp.seekg(0, std::ios::end);
    const std::streampos file_size = fp.tellg();
    fp.seekg(0, std::ios::beg);

    // Allocate memory
    std::string file_content;
    file_content.reserve(file_size);

    // Read the file
    file_content.assign((std::istreambuf_iterator<char>(fp)),
            std::istreambuf_iterator<char>());
    fp.close();
    ryml::EventHandlerTree event_handler = {};
    ryml::Parser parser(&event_handler, ryml::ParserOptions().locations(true));
    const ryml::Tree top_level_yaml =
            ryml::parse_in_place(&parser, file_content.data());
    INFO("Loading network YAML information from file.\n");
    // NOLINTNEXT(misc-include-cleaner)
    ryml::Tree tree = ryml::parse_in_place(file_content.data());
    INFO("Network YAML information loaded from file.\n");

    const ryml::ConstNodeRef yaml_node = tree.rootref();
    if (yaml_node.is_map())
    {
        if (yaml_node.find_child("network").invalid())
        {
            throw DescriptionParsingError(
                    "No top-level 'network' section defined", parser,
                    yaml_node);
        }
        SpikingNetwork net = description_parse_network_section_yaml(
                parser, yaml_node["network"]);
        if (yaml_node.find_child("mappings").invalid())
        {
            throw DescriptionParsingError(
                    "No 'mappings' section defined", parser, yaml_node);
        }
        description_parse_mapping_section_yaml(
                parser, yaml_node["mappings"], arch, net);

        return net;
    }
    throw DescriptionParsingError(
            "Mapped network file has invalid format", parser, yaml_node);
}

sanafe::SpikingNetwork sanafe::description_parse_network_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef net_node)
{
    std::string net_name;
    if (!net_node.find_child("name").invalid())
    {
        net_node["name"] >> net_name;
        if (net_name.find('[') != std::string::npos)
        {
            throw DescriptionParsingError(
                    "Multiple networks not supported", parser, net_node);
        }
    }
    else
    {
        INFO("Warning: No network name given; leaving name empty.\n");
    }

    INFO("Parsing network: %s\n", net_name.c_str());

    SpikingNetwork new_net(std::move(net_name));
    description_parse_neuron_group_section_yaml(
            parser, net_node["groups"], new_net);
    description_parse_edges_section_yaml(parser, net_node["edges"], new_net);

    return new_net;
}

void sanafe::description_parse_neuron_group_section_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef groups_node,
        SpikingNetwork &net)
{
    INFO("Parsing neuron groups.\n");
    if (groups_node.is_seq())
    {
        for (const auto &group : groups_node)
        {
            description_parse_group(parser, group, net);
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Neuron group section does not define a list of groups", parser,
                groups_node);
    }
}

void sanafe::description_parse_edges_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef edges_node, SpikingNetwork &net)
{
    INFO("Parsing edges section.\n");
    if (edges_node.is_seq())
    {
        for (const auto list_entry : edges_node)
        {
            for (const auto edge_node : list_entry)
            {
                std::string edge_description;
                // NOLINTNEXT(misc-include-cleaner)
                edge_node >> ryml::key(edge_description);
                description_parse_edge(
                        edge_description, parser, edge_node, net);
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Edges section does not define a list of edges", parser,
                edges_node);
    }
}

void sanafe::description_parse_group(const ryml::Parser &parser,
        const ryml::ConstNodeRef neuron_group_node, SpikingNetwork &net)
{
    auto group_name = description_required_field<std::string>(
            parser, neuron_group_node, "name");
    INFO("Parsing neuron group: %s\n", group_name.c_str());

    if (neuron_group_node.find_child("neurons").invalid())
    {
        throw DescriptionParsingError(
                "No neurons section defined.", parser, neuron_group_node);
    }
    const auto &neurons_node = neuron_group_node["neurons"];
    const size_t neuron_count = description_count_neurons(parser, neurons_node);

    NeuronConfiguration default_neuron_config{};
    if (!neuron_group_node.find_child("attributes").invalid())
    {
        TRACE1(DESCRIPTION, "Parsing neuron group attributes\n");
        default_neuron_config = description_parse_neuron_attributes_yaml(
                parser, neuron_group_node["attributes"]);
    }
    NeuronGroup &group = net.create_neuron_group(
            std::move(group_name), neuron_count, default_neuron_config);
    TRACE1(DESCRIPTION, "Parsing neuron section\n");
    description_parse_neuron_section_yaml(parser, neurons_node, group);
}

size_t sanafe::description_count_neurons(
        const ryml::Parser &parser, const ryml::ConstNodeRef neuron_node)
{
    size_t neuron_count{0UL};

    if (neuron_node.is_seq())
    {
        for (const auto &neuron_entry : neuron_node)
        {
            if (neuron_entry.is_map() || neuron_entry.is_seq())
            {
                for (const auto neuron_description : neuron_entry)
                {
                    std::string id;
                    neuron_description >> ryml::key(id);
                    const bool is_range = (id.find("..") != std::string::npos);
                    if (is_range)
                    {
                        const auto range = description_parse_range_yaml(id);
                        neuron_count += (range.second - range.first) + 1;
                    }
                    else // if entry defines a single neuron
                    {
                        ++neuron_count;
                    }
                }
            }
            else
            {
                std::string id;
                neuron_entry >> id;
                const bool is_range = (id.find("..") != std::string::npos);
                if (is_range)
                {
                    const auto range = description_parse_range_yaml(id);
                    neuron_count += (range.second - range.first) + 1;
                }
                else // if entry defines a single neuron
                {
                    ++neuron_count;
                }
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format, should be list", parser, neuron_node);
    }

    INFO("Counted %zu neurons\n", neuron_count);
    return neuron_count;
}

void sanafe::description_parse_neuron_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef neuron_node, NeuronGroup &neuron_group)
{
    if (neuron_node.is_seq())
    {
        for (const auto &list_entry : neuron_node)
        {
            for (const auto &neuron_description : list_entry)
            {
                // Iterate, but there should only be one mapping per list entry
                std::string id;
                neuron_description >> ryml::key(id);
                description_parse_neuron(
                        id, parser, neuron_description, neuron_group);
            }
        }
    }
    else
    {
        throw DescriptionParsingError(
                "Invalid neuron format, should be list", parser, neuron_node);
    }
}

void sanafe::description_parse_neuron(const std::string &id,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        NeuronGroup &neuron_group)
{
    std::pair<size_t, size_t> range;
    TRACE1(DESCRIPTION, "Parsing neuron(s): %s\n", id.c_str());
    const NeuronConfiguration config = description_parse_neuron_attributes_yaml(
            parser, attributes, neuron_group.default_neuron_config);
    const bool is_range = (id.find("..") != std::string::npos);
    if (is_range)
    {
        range = description_parse_range_yaml(id);
        for (size_t instance = range.first; instance <= range.second;
                ++instance)
        {
            Neuron &n = neuron_group.neurons.at(instance);
            n.set_attributes(config);
        }
    }
    else
    {
        const size_t nid = std::stoull(id);
        Neuron &n = neuron_group.neurons.at(nid);
        n.set_attributes(config);
    }
}

sanafe::NeuronConfiguration sanafe::description_parse_neuron_attributes_yaml(
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes,
        const NeuronConfiguration &default_template)
{
    NeuronConfiguration neuron_template = default_template;

    if (attributes.is_seq())
    {
        // Ordered list format, recursively parse attributes for each element
        for (const auto &attribute : attributes)
        {
            neuron_template = description_parse_neuron_attributes_yaml(
                    parser, attribute, neuron_template);
        }
        return neuron_template;
    }

    if (!attributes.find_child("log_potential").invalid())
    {
        bool log_potential = false;
        attributes["log_potential"] >> log_potential;
        neuron_template.log_potential = log_potential;
    }
    if (!attributes.find_child("log_spikes").invalid())
    {
        bool log_spikes = false;
        attributes["log_spikes"] >> log_spikes;
        neuron_template.log_spikes = log_spikes;
    }
    if (!attributes.find_child("force_synapse_update").invalid())
    {
        bool force_synapse_update = false;
        attributes["force_synapse_update"] >> force_synapse_update;
        neuron_template.force_synapse_update = force_synapse_update;
    }
    if (!attributes.find_child("force_dendrite_update").invalid())
    {
        bool force_dendrite_update = false;
        attributes["force_dendrite_update"] >> force_dendrite_update;
        neuron_template.force_dendrite_update = force_dendrite_update;
    }
    if (!attributes.find_child("force_soma_update").invalid())
    {
        bool force_soma_update = false;
        attributes["force_soma_update"] >> force_soma_update;
        neuron_template.force_soma_update = force_soma_update;
    }
    if (!attributes.find_child("synapse_hw_name").invalid())
    {
        std::string synapse_hw_name;
        attributes["synapse_hw_name"] >> synapse_hw_name;
        neuron_template.default_synapse_hw_name = synapse_hw_name;
    }
    if (!attributes.find_child("dendrite_hw_name").invalid())
    {
        std::string dendrite_hw_name;
        attributes["dendrite_hw_name"] >> dendrite_hw_name;
        neuron_template.dendrite_hw_name = dendrite_hw_name;
    }
    if (!attributes.find_child("soma_hw_name").invalid())
    {
        std::string soma_hw_name;
        attributes["soma_hw_name"] >> soma_hw_name;
        neuron_template.soma_hw_name = soma_hw_name;
    }

    // Parse and add shared attributes, which are defined alongside attributes
    auto model_attributes =
            description_parse_model_attributes_yaml(parser, attributes);
    for (auto &[key, attribute] : model_attributes)
    {
        attribute.forward_to_dendrite = true;
        attribute.forward_to_soma = true;
        neuron_template.model_attributes[key] = attribute;
    }
    // Parse h/w unit specific model attributes
    if (!attributes.find_child("dendrite").invalid())
    {
        auto dendrite_attributes = description_parse_model_attributes_yaml(
                parser, attributes["dendrite"]);
        for (auto &[key, attribute] : dendrite_attributes)
        {
            attribute.forward_to_synapse = false;
            attribute.forward_to_soma = false;
            neuron_template.model_attributes[key] = attribute;
        }
    }
    if (!attributes.find_child("soma").invalid())
    {
        auto soma_attributes = description_parse_model_attributes_yaml(
                parser, attributes["soma"]);
        for (auto &[key, attribute] : soma_attributes)
        {
            attribute.forward_to_synapse = false;
            attribute.forward_to_dendrite = false;
            neuron_template.model_attributes[key] = attribute;
        }
    }

    return neuron_template;
}

static std::string_view description_trim_whitespace(
        const std::string_view input)
{
    constexpr auto whitespace = " \t\n\r";
    auto start = input.find_first_not_of(whitespace);
    auto end = input.find_last_not_of(whitespace);
    if ((start == std::string::npos) || (end == std::string::npos))
    {
        return "";
    }
    return input.substr(start, end - start + 1);
}

std::tuple<sanafe::NeuronAddress, sanafe::NeuronAddress>
sanafe::description_parse_edge_description(const std::string_view &description)
{
    auto arrow_pos = description.find("->");
    if (arrow_pos == std::string::npos)
    {
        throw std::runtime_error(
                "Edge is not formatted correctly: " + std::string(description));
    }

    const std::string_view source_part =
            description_trim_whitespace(description.substr(0, arrow_pos));
    const std::string_view target_part =
            description_trim_whitespace(description.substr(arrow_pos + 2));

    const auto source_dot_pos = source_part.find('.');
    const auto target_dot_pos = target_part.find('.');

    const bool source_neuron_defined = (source_dot_pos != std::string::npos);
    const bool target_neuron_defined = (target_dot_pos != std::string::npos);
    if (source_neuron_defined && !target_neuron_defined)
    {
        throw std::runtime_error(
                "No target neuron defined in edge:" + std::string(description));
    }
    if (target_neuron_defined && !source_neuron_defined)
    {
        throw std::runtime_error(
                "No target neuron defined in edge:" + std::string(description));
    }

    NeuronAddress source;
    source.group_name = source_part.substr(0, source_dot_pos);
    if (source_neuron_defined)
    {
        source.neuron_offset = std::stoull(
                std::string(source_part.substr(source_dot_pos + 1)));
    }

    NeuronAddress target;
    target.group_name = target_part.substr(0, target_dot_pos);
    if (target_neuron_defined)
    {
        target.neuron_offset = std::stoull(
                std::string(target_part.substr(target_dot_pos + 1)));
    }

    return std::make_tuple(source, target);
}

void sanafe::description_parse_edge(const std::string &description,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node,
        SpikingNetwork &net)
{
    // Description has format src_group.src_neuron -> tgt_group.tgt_neuron
    const auto [source_address, target_address] =
            description_parse_edge_description(description);

    const bool is_hyper_edge = !source_address.neuron_offset.has_value();
    if (is_hyper_edge)
    {
        description_parse_hyperedge(
                source_address, target_address, parser, attributes_node, net);
    }
    else
    {
        description_parse_neuron_connection(
                source_address, target_address, parser, attributes_node, net);
    }
}

void sanafe::description_parse_neuron_connection(
        const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        const ryml::ConstNodeRef attributes_node, SpikingNetwork &net)
{
    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &source_group = net.groups.at(source_address.group_name);
    if (source_address.neuron_offset >= source_group.neurons.size())
    {
        std::string error =
                "Invalid source neuron id: " + source_address.group_name;
        if (source_address.neuron_offset.has_value())
        {
            error += "." + std::to_string(source_address.neuron_offset.value());
        }
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &source_neuron =
            source_group.neurons[source_address.neuron_offset.value()];

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    NeuronGroup &target_group = net.groups.at(target_address.group_name);

    if (target_address.neuron_offset >= target_group.neurons.size())
    {
        const std::string error =
                "Invalid target neuron id: " + target_address.group_name + "." +
                std::to_string(target_address.neuron_offset.value());
        throw DescriptionParsingError(error, parser, attributes_node);
    }
    Neuron &target_neuron =
            target_group.neurons.at(target_address.neuron_offset.value());

    const size_t edge_idx = source_neuron.connect_to_neuron(target_neuron);
    Connection &edge = source_neuron.edges_out[edge_idx];
    description_parse_edge_attributes(edge, parser, attributes_node);
}

void sanafe::description_parse_hyperedge(const NeuronAddress &source_address,
        const NeuronAddress &target_address, const ryml::Parser &parser,
        const ryml::ConstNodeRef hyperedge_node, SpikingNetwork &net)
{
    if (net.groups.find(source_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid source neuron group:" + source_address.group_name;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }
    NeuronGroup &source_group = net.groups.at(source_address.group_name);

    if (net.groups.find(target_address.group_name) == net.groups.end())
    {
        const std::string error =
                "Invalid target neuron group:" + target_address.group_name;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }

    NeuronGroup &target_group = net.groups.at(target_address.group_name);

    // First parse the lists of attributes, this is common for all hyperedge
    //  connectivity
    const auto type = description_required_field<std::string>(
            parser, hyperedge_node, "type");

    std::map<std::string, std::vector<ModelAttribute>> attribute_lists{};
    if (type == "conv2d")
    {
        Conv2DParameters convolution{};
        for (const auto &attribute : hyperedge_node.children())
        {
            if (attribute.key() == "input_height")
            {
                attribute >> convolution.input_height;
            }
            else if (attribute.key() == "input_width")
            {
                attribute >> convolution.input_width;
            }
            else if (attribute.key() == "input_channels")
            {
                attribute >> convolution.input_channels;
            }
            else if (attribute.key() == "kernel_width")
            {
                attribute >> convolution.kernel_width;
            }
            else if (attribute.key() == "kernel_height")
            {
                attribute >> convolution.kernel_height;
            }
            else if (attribute.key() == "kernel_count")
            {
                attribute >> convolution.kernel_count;
            }
            else if (attribute.key() == "stride_width")
            {
                attribute >> convolution.stride_width;
            }
            else if (attribute.key() == "stride_height")
            {
                attribute >> convolution.stride_height;
            }
            else if (attribute.key() == "synapse")
            {
                for (const auto &synapse_attribute_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelAttribute> attribute_list;
                    for (const auto &model_attribute_node :
                            synapse_attribute_node)
                    {
                        ModelAttribute value = description_parse_attribute_yaml(
                                parser, model_attribute_node);
                        value.forward_to_dendrite = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(std::move(value));
                    }
                    std::string attribute_name;
                    synapse_attribute_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = std::move(attribute_list);
                }
            }
            else if (attribute.key() == "dendrite")
            {
                for (const auto &dendrite_attribute_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelAttribute> attribute_list;
                    for (const auto &model_attribute_node :
                            dendrite_attribute_node)
                    {
                        ModelAttribute value = description_parse_attribute_yaml(
                                parser, model_attribute_node);
                        value.forward_to_synapse = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(std::move(value));
                    }
                    std::string attribute_name;
                    dendrite_attribute_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = std::move(attribute_list);
                }
            }
            else if (attribute.key() != "type")
            {
                // TODO: refactor
                std::vector<ModelAttribute> attribute_list;
                for (const auto &model_attribute_node : attribute)
                {
                    ModelAttribute value = description_parse_attribute_yaml(
                            parser, model_attribute_node);
                    attribute_list.push_back(std::move(value));
                }
                std::string attribute_name;
                attribute >> ryml::key(attribute_name);
                attribute_lists[attribute_name] = std::move(attribute_list);
            }
        }

        source_group.connect_neurons_conv2d(
                target_group, attribute_lists, convolution);
    }
    else if (type == "dense")
    {
        for (const auto &attribute : hyperedge_node.children())
        {
            // TODO: refactor
            if (attribute.key() == "synapse")
            {
                for (const auto &synapse_attribute_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelAttribute> attribute_list;
                    for (const auto &model_attribute_node :
                            synapse_attribute_node)
                    {
                        ModelAttribute value = description_parse_attribute_yaml(
                                parser, model_attribute_node);
                        value.forward_to_dendrite = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(std::move(value));
                    }
                    std::string attribute_name;
                    synapse_attribute_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = std::move(attribute_list);
                }
            }
            else if (attribute.key() == "dendrite")
            {
                for (const auto &dendrite_attribute_node : attribute)
                {
                    // TODO: refactor
                    std::vector<ModelAttribute> attribute_list;
                    for (const auto &model_attribute_node :
                            dendrite_attribute_node)
                    {
                        ModelAttribute value = description_parse_attribute_yaml(
                                parser, model_attribute_node);
                        value.forward_to_synapse = false;
                        value.forward_to_soma = false;
                        attribute_list.push_back(std::move(value));
                    }
                    std::string attribute_name;
                    dendrite_attribute_node >> ryml::key(attribute_name);
                    attribute_lists[attribute_name] = std::move(attribute_list);
                }
            }
            else if (attribute.key() != "type")
            {
                // TODO: refactor
                std::vector<ModelAttribute> attribute_list;
                for (const auto &model_attribute_node : attribute)
                {
                    ModelAttribute value = description_parse_attribute_yaml(
                            parser, model_attribute_node);
                    attribute_list.push_back(std::move(value));
                }
                std::string attribute_name;
                attribute >> ryml::key(attribute_name);
                attribute_lists[attribute_name] = std::move(attribute_list);
            }
        }
        source_group.connect_neurons_dense(target_group, attribute_lists);
    }
    else if (type == "sparse")
    {
    }
    else
    {
        const std::string error = "Invalid hyperedge type: " + type;
        throw DescriptionParsingError(error, parser, hyperedge_node);
    }
}

void sanafe::description_parse_edge_attributes(Connection &edge,
        const ryml::Parser &parser, const ryml::ConstNodeRef attributes_node)
{
    if (attributes_node.is_seq())
    {
        for (const auto &attribute : attributes_node)
        {
            description_parse_edge_attributes(edge, parser, attribute);
        }

        return;
    }

    if (!attributes_node.find_child("synapse").invalid())
    {
        auto synapse_attributes = description_parse_model_attributes_yaml(
                parser, attributes_node["synapse"]);
        for (auto &[key, attribute] : synapse_attributes)
        {
            attribute.forward_to_dendrite = false;
            attribute.forward_to_soma = false;
            edge.synapse_attributes[key] = attribute;
        }
    }
    if (!attributes_node.find_child("dendrite").invalid())
    {
        auto dendrite_attributes = description_parse_model_attributes_yaml(
                parser, attributes_node["dendrite"]);
        for (auto &[key, attribute] : dendrite_attributes)
        {
            attribute.forward_to_synapse = false;
            attribute.forward_to_soma = false;
            edge.dendrite_attributes[key] = attribute;
        }
    }

    const auto shared_model_attributes =
            description_parse_model_attributes_yaml(parser, attributes_node);
    for (const auto &[key, attribute] : shared_model_attributes)
    {
        if ((key != "synapse") && (key != "dendrite") && (key != "soma"))
        {
            TRACE2(DESCRIPTION, "Adding con attribute:%s\n", key.c_str());
            edge.synapse_attributes[key] = attribute;
            edge.dendrite_attributes[key] = attribute;
        }
    }
}

void sanafe::description_parse_mapping_section_yaml(const ryml::Parser &parser,
        const ryml::ConstNodeRef mappings_node, Architecture &arch,
        SpikingNetwork &net)
{
    INFO("Parsing mapping section.\n");
    if (!mappings_node.is_seq())
    {
        throw DescriptionParsingError(
                "Mappings must be given as a sequence / list. Each list entry "
                "must first use the neuron address as a key, followed by the "
                "core address.\nE.g.: "
                "- G.n: {core: 1.1}\nThis maps group G's neuron 'n' "
                "to Tile 1's Core 1, specifically using soma unit 'foo'.)",
                parser, mappings_node);
    }

    for (const auto &mapping : mappings_node)
    {
        if (!mapping.is_map())
        {
            throw DescriptionParsingError(
                    "Expected mapping to be defined in the format: "
                    "<group>.<neuron>: [<attributes>]",
                    parser, mapping);
        }

        int entries = 0;
        // Ordered mappings, so each list entry contains a mapping with a single
        //  key. Only look up the first entry.
        for (const auto &mapping_entry : mapping)
        {
            description_parse_mapping(parser, mapping_entry, arch, net);
            ++entries;
            if (entries > 1)
            {
                throw DescriptionParsingError("Should be one entry per mapping",
                        parser, mapping_entry);
            }
        }
    }
}

void sanafe::description_parse_mapping(const ryml::Parser &parser,
        const ryml::ConstNodeRef mapping_info, Architecture &arch,
        SpikingNetwork &net)
{
    std::string neuron_address;
    mapping_info >> ryml::key(neuron_address);
    const auto dot_pos = neuron_address.find_first_of('.');
    const bool neuron_defined = dot_pos != std::string::npos;

    const std::string group_name = neuron_address.substr(0, dot_pos);
    NeuronGroup &group = net.groups.at(group_name);

    size_t start_id = 0;
    size_t end_id = 0;
    if (neuron_defined)
    {
        const std::string neuron_str = neuron_address.substr(dot_pos + 1);
        if (net.groups.find(group_name) == net.groups.end())
        {
            const std::string error = "Invalid neuron group:" + group_name;
            throw DescriptionParsingError(error, parser, mapping_info);
        }

        if (neuron_str.find("..") != std::string::npos)
        {
            std::tie(start_id, end_id) =
                    description_parse_range_yaml(neuron_str);
        }
        else
        {
            start_id = std::stoull(neuron_str);
            end_id = start_id;
        }
    }
    else
    {
        // No neuron given so map all neurons in the group
        start_id = 0UL;
        assert(group.neurons.size() > 0);
        end_id = group.neurons.size() - 1UL;
    }

    //INFO("Mapping neuron: %s.%zu\n", group_name.c_str(), neuron_id);
    for (size_t neuron_offset = start_id; neuron_offset <= end_id;
            ++neuron_offset)
    {
        if (neuron_offset >= group.neurons.size())
        {
            std::string error = "Invalid neuron id: ";
            error += group_name;
            error += '.';
            error += std::to_string(neuron_offset);
            throw DescriptionParsingError(error, parser, mapping_info);
        }

        // Get any mapping attributes or configuration
        std::string core_address;
        Neuron &n = group.neurons[neuron_offset];
        description_parse_mapping_info(parser, mapping_info, n, core_address);

        // Get pointers to the h/w we're mapping to
        const auto dot_pos = core_address.find('.');
        const size_t tile_id = std::stoull(core_address.substr(0, dot_pos));
        const size_t core_offset_within_tile =
                std::stoull(core_address.substr(dot_pos + 1));

        if (tile_id >= arch.tiles.size())
        {
            throw DescriptionParsingError(
                    "Tile ID >= tile count", parser, mapping_info);
        }
        TileConfiguration &tile = arch.tiles[tile_id];
        if (core_offset_within_tile >= tile.cores.size())
        {
            throw DescriptionParsingError(
                    "Core ID >= core count", parser, mapping_info);
        }
        const CoreConfiguration &core = tile.cores[core_offset_within_tile];
        n.map_to_core(core);
    }
}

void sanafe::description_parse_mapping_info(const ryml::Parser &parser,
        const ryml::ConstNodeRef info, Neuron &n, std::string &core_name)
{
    if (info.is_seq())
    {
        for (const auto &field : info)
        {
            description_parse_mapping_info(parser, field, n, core_name);
        }
    }
    else if (!info.is_map())
    {
        throw DescriptionParsingError(
                "Expected attributes to be map", parser, info);
    }
    else
    {
        if (!info.find_child("synapse").invalid())
        {
            info["synapse"] >> n.default_synapse_hw_name;
        }
        if (!info.find_child("dendrite").invalid())
        {
            info["dendrite"] >> n.dendrite_hw_name;
        }
        if (!info.find_child("soma").invalid())
        {
            info["soma"] >> n.soma_hw_name;
        }
        if (!info.find_child("core").invalid())
        {
            info["core"] >> core_name;
        }
    }
}

void sanafe::description_write_network_yaml(
        const std::filesystem::path path, const sanafe::SpikingNetwork &network)
{
    std::ifstream previous_content(path);

    // Create a new YAML tree
    ryml::Tree tree;
    ryml::NodeRef root = tree.rootref();

    // Try to read existing content if file exists and is not empty
    const bool file_empty =
            (previous_content.peek() == std::ifstream::traits_type::eof());

    if (!previous_content.is_open() || file_empty)
    {
        // Read existing YAML content
        const std::string existing_content(
                (std::istreambuf_iterator<char>(previous_content)),
                std::istreambuf_iterator<char>());

        // Parse existing content
        tree = ryml::parse_in_arena(existing_content.c_str());
        root = tree.rootref();
        previous_content.close();

        // Remove the existing network info if it exists
        if (root.has_child("network"))
        {
            root.remove_child("network");
        }
    }
    else
    {
        // Initialize with empty document
        root |= ryml::MAP;
    }

    std::ofstream fp(path);
    if (!fp.is_open())
    {
        throw std::system_error(std::make_error_code(std::errc::io_error),
                "Failed to open YAML file for writing: " + path.string());
    }
    // Note we need to keep all the generated strings alive as long as we are
    //  dealing with this YAML file. RapidYAML only uses views of strings for
    //  speed, and doesn't copy the string contents
    std::list<std::string> strings{};
    // Add network section
    description_serialize_network_to_yaml(root, network, strings);

    // Convert to string and write to file
    std::ostringstream ss;
    ss << tree;
    fp << ss.str();
    fp.close();
}

ryml::NodeRef sanafe::description_serialize_network_to_yaml(ryml::NodeRef root,
        const sanafe::SpikingNetwork &network, std::list<std::string> &strings)
{
    auto network_node = root["network"];
    network_node |= ryml::MAP;
    if (network.name.empty())
    {
        network_node["name"] << " ";
    }
    else
    {
        network_node["name"] << network.name;
    }

    // Add neuron groups
    auto groups_node = network_node["groups"];
    groups_node |= ryml::SEQ; // NOLINT(misc-include-cleaner)

    for (const auto &[name, group] : network.groups)
    {
        description_serialize_neuron_group_to_yaml(groups_node, group, strings);
    }

    // Add edges (connections)
    auto edges_node = network_node["edges"];
    edges_node |= ryml::SEQ; // NOLINT(misc-include-cleaner)

    // Iterate through all neurons and their connections
    for (const auto &[group_name, group] : network.groups)
    {
        for (const auto &neuron : group.neurons)
        {
            for (const auto &connection : neuron.edges_out)
            {
                ryml::NodeRef edge_map = edges_node.append_child();
                edge_map |= ryml::MAP; // NOLINT(misc-include-cleaner)

                // Create edge description (source -> destination)
                const std::string edge_description =
                        connection.pre_neuron.group_name + "." +
                        std::to_string(
                                connection.pre_neuron.neuron_offset.value()) +
                        " -> " + connection.post_neuron.group_name + "." +
                        std::to_string(
                                connection.post_neuron.neuron_offset.value());

                const std::string &ref = strings.emplace_back(edge_description);
                ryml::NodeRef edge_node = edge_map[ref.c_str()];
                edge_node |= ryml::MAP; // NOLINT(misc-include-cleaner)
                // For conciseness use flow style outputs for edge attributes
                edge_node |= ryml::FLOW_SL; // NOLINT(misc-include-cleaner)
                // For now assume there are no default connection attributes
                const std::map<std::string, ModelAttribute> default_attributes{};
                description_serialize_model_attributes_to_yaml(edge_node,
                        connection.synapse_attributes, default_attributes);

                // TODO: support synapse-specific attributes
                // if (!connection.synapse_attributes.empty())
                // {
                //     auto synapse_node = edge_node["synapse"];
                //     synapse_node |= ryml::MAP;
                //     description_serialize_model_attributes_to_yaml(
                //             synapse_node, connection.synapse_attributes);
                // }
                // // Add dendrite-specific attributes
                // if (!connection.dendrite_attributes.empty())
                // {
                //     auto dendrite_node = edge_node["dendrite"];
                //     dendrite_node |= ryml::MAP;
                //     description_serialize_model_attributes_to_yaml(
                //             dendrite_node, connection.dendrite_attributes);
                // }
            }
        }
    }

    return network_node;
}

ryml::NodeRef sanafe::description_serialize_neuron_group_to_yaml(
        ryml::NodeRef parent, const sanafe::NeuronGroup &group,
        std::list<std::string> &strings)
{
    auto group_node = parent.append_child();
    group_node |= ryml::MAP;
    group_node["name"] << group.name;

    // Add attributes if they exist
    auto attr_node = group_node["attributes"];
    attr_node |= ryml::MAP;

    // Add model attributes if they exist
    const std::map<std::string, ModelAttribute> no_default_attributes{};
    if (!group.default_neuron_config.model_attributes.empty())
    {
        description_serialize_model_attributes_to_yaml(attr_node,
                group.default_neuron_config.model_attributes,
                no_default_attributes);
    }

    // Add neurons
    auto neurons_node = group_node["neurons"];
    neurons_node |= ryml::SEQ;

    // Group neurons with identical configurations to reduce output size
    std::vector<std::tuple<size_t, size_t>> neuron_runs;
    size_t run_start{0UL};
    auto prev_neuron = group.neurons.begin();
    for (auto neuron_it = group.neurons.begin();
            neuron_it != group.neurons.end(); ++neuron_it)
    {
        if (neuron_it->model_attributes != prev_neuron->model_attributes)
        {
            neuron_runs.emplace_back(run_start, prev_neuron->offset);
            INFO("adding new run %zu..%zu\n", run_start, prev_neuron->offset);
            // Set up the next run of unique neurons
            run_start = neuron_it->offset;
        }
        prev_neuron = neuron_it;
    }
    neuron_runs.emplace_back(run_start, prev_neuron->offset);
    INFO("adding new run %zu..%zu\n", run_start, prev_neuron->offset);

    for (const auto &neuron_run : neuron_runs)
    {
        auto [start_offset, end_offset] = neuron_run;

        auto neuron_map = neurons_node.append_child();
        neuron_map |= ryml::MAP; // NOLINT(misc-include-cleaner)

        const Neuron &neuron = group.neurons[start_offset];
        std::string neuron_description = std::to_string(start_offset);
        if (end_offset != start_offset)
        {
            neuron_description += ".." + std::to_string(end_offset);
        }
        const std::string &ref = strings.emplace_back(neuron_description);
        auto neuron_node = neuron_map[ref.c_str()];

        // Add model attributes if they exist and differ from group defaults
        neuron_node |= ryml::MAP; // NOLINT(misc-include-cleaner)
        neuron_node |= ryml::FLOW_SL; // NOLINT(misc-include-cleaner)
        if (!neuron.model_attributes.empty())
        {
            description_serialize_model_attributes_to_yaml(neuron_node,
                    neuron.model_attributes,
                    group.default_neuron_config.model_attributes);
        }
    }

    return group_node;
}

ryml::NodeRef sanafe::description_serialize_model_attributes_to_yaml(
        ryml::NodeRef parent,
        const std::map<std::string, sanafe::ModelAttribute> &attributes,
        const std::map<std::string, sanafe::ModelAttribute> &default_values)
{
    // Add all attributes to the parent node
    for (const auto &[key, attribute] : attributes)
    {
        TRACE1(DESCRIPTION, "Adding attribute %s\n", key.c_str());
        // TODO: support model-specific attributes
        if (key == "synapse" || key == "dendrite" || key == "soma")
        {
            continue;
        }

        if ((default_values.find(key) == default_values.end()) ||
                (default_values.at(key) != attribute))
        {
            // Its safe to reference into these strings; they will still valid
            //  over the duration of the save
            description_serialize_variant_value_to_yaml(
                    parent[key.c_str()], attribute.value);
        }
    }

    return parent;
}

ryml::NodeRef sanafe::description_serialize_variant_value_to_yaml(
        ryml::NodeRef node,
        const std::variant<bool, int, double, std::string,
                std::vector<sanafe::ModelAttribute>> &value)
{
    std::visit(
            [&node](auto &&arg) {
                using T = std::decay_t<decltype(arg)>;

                if constexpr (std::is_same_v<T, int>)
                {
                    node << arg;
                }
                else if constexpr (std::is_same_v<T, double>)
                {
                    node << arg;
                }
                else if constexpr (std::is_same_v<T, bool>)
                {
                    node << arg;
                }
                else if constexpr (std::is_same_v<T, std::string>)
                {
                    node << arg;
                }
                else if constexpr (std::is_same_v<T,
                                           std::vector<ModelAttribute>>)
                {
                    // Handle list of attributes
                    node |= ryml::SEQ;

                    for (const ModelAttribute &attribute : arg)
                    {
                        auto child = node.append_child();

                        if (!attribute.name.has_value() ||
                                attribute.name.value().empty())
                        {
                            // Unnamed attribute - directly serialize its value
                            description_serialize_variant_value_to_yaml(
                                    child, attribute.value);
                        }
                        else
                        {
                            // Named attribute - create a map with the name as
                            //  key
                            child |= ryml::MAP; // NOLINT(misc-include-cleaner)
                            auto attribute_node =
                                    child[attribute.name.value().c_str()];
                            description_serialize_variant_value_to_yaml(
                                    attribute_node, attribute.value);
                        }
                    }
                }
            },
            value);

    return node;
}

void sanafe::description_write_mappings_yaml(
        const std::filesystem::path path, const sanafe::SpikingNetwork &network)
{
    std::ifstream previous_content(path);

    // Create a new YAML tree
    ryml::Tree tree;
    ryml::NodeRef root = tree.rootref();

    // Try to read existing content if file exists and is not empty
    const bool file_empty =
            (previous_content.peek() == std::ifstream::traits_type::eof());
    if (!previous_content.is_open() || !file_empty)
    {
        // Read existing YAML content
        const std::string existing_content(
                (std::istreambuf_iterator<char>(previous_content)),
                std::istreambuf_iterator<char>());

        // Parse existing content
        tree = ryml::parse_in_arena(existing_content.c_str());
        root = tree.rootref();
        previous_content.close();

        // Remove the existing network info if it exists
        if (root.has_child("mappings"))
        {
            root.remove_child("mappings");
        }
    }
    else
    {
        // Initialize with empty document
        root |= ryml::MAP; // NOLINT(misc-include-cleaner)
    }

    std::ofstream fp(path);
    if (!fp.is_open())
    {
        throw std::system_error(std::make_error_code(std::errc::io_error),
                "Failed to open YAML file for writing: " + path.string());
    }
    std::list<std::string> strings{};

    // Add mappings section
    auto mappings_node = root["mappings"];
    mappings_node |= ryml::SEQ;

    // Save all mappings
    std::vector<std::reference_wrapper<const Neuron>> all_neurons;

    // Collect all neurons from all groups
    for (const auto &group : network.groups)
    {
        const auto &neurons = group.second.neurons;
        all_neurons.insert(all_neurons.end(), neurons.begin(), neurons.end());
    }

    // Sort by mapping_order
    std::sort(all_neurons.begin(), all_neurons.end(),
            [](const Neuron &a, const Neuron &b) {
                return a.mapping_order < b.mapping_order;
            });

    // Now write mappings
    for (const Neuron &neuron : all_neurons)
    {
        if (!neuron.core_address.has_value())
        {
            INFO("Error: Neuron (nid:%s.%zu) not mapped, can't save.\n",
                    neuron.parent_group_name.c_str(), neuron.offset);
            throw std::runtime_error("Error: Neuron not mapped, can't save.");
        }

        auto mapping_entry = mappings_node.append_child();
        mapping_entry |= ryml::MAP; // NOLINT(misc-include-cleaner)

        std::string neuron_addr;
        neuron_addr =
                neuron.parent_group_name + "." + std::to_string(neuron.offset);

        const std::string &neuron_ref = strings.emplace_back(neuron_addr);
        auto mapping_info = mapping_entry[neuron_ref.c_str()];
        mapping_info |= ryml::MAP; // NOLINT(misc-include-cleaner)

        // Add core address
        const std::string core_address =
                std::to_string(neuron.core_address->parent_tile_id) + "." +
                std::to_string(neuron.core_address->id);
        const std::string &core_ref = strings.emplace_back(core_address);
        mapping_info["core"] << core_ref;

        // Add hardware component names if consistently used by all neurons in
        //  range
        if (!neuron.soma_hw_name.empty())
        {
            mapping_info["soma"] << neuron.soma_hw_name;
        }

        if (!neuron.dendrite_hw_name.empty())
        {
            mapping_info["dendrite"] << neuron.dendrite_hw_name;
        }

        if (!neuron.default_synapse_hw_name.empty())
        {
            mapping_info["synapse"] << neuron.default_synapse_hw_name;
        }
    }

    // Convert to string and write to file
    std::stringstream ss;
    ss << tree;
    fp << ss.str();
    fp.close();
}
