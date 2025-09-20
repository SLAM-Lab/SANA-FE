# Capstone code for live demo
timesteps = 0
if run_alive:
    while True:
        if timesteps > 0:
        sim.run(timesteps, heartbeat=100)
        sim.get_run_summary()
        print("Enter timesteps to run: ", end="")
        user_in = input()

        if user_in == "q" or user_in == "quit":
            break
        if user_in.startswith("u"):
            try:
                group_id = int(user_in[2])
            except ValueError:
                print(f"Error: Expected int. Got \"{user_in[2]}\".")
                exit(1)

            try:
                n_id = int(user_in[4])
            except ValueError:
                print(f"Error: Expected int. Got \"{user_in[4]}\".")
                exit(1)

            user_in = user_in[6:]
            kwargs = user_in.split(" ")
            # print(group_id, n_id, kwargs, len(kwargs))
            #sim.update_neuron(group_id, n_id, kwargs, len(kwargs))

            timesteps = 0
            continue
        if user_in.startswith("s"):
            try:
                group_id = int(user_in[2])
            except ValueError:
                print(f"Error: Expected int. Got \"{user_in[2]}\".")
                exit(1)

            #print(sim.get_status(group_id))

            timesteps = 0
            continue

        try:
            timesteps = int(user_in)
        except ValueError:
            print(f"Error: Expected int. Got {user_in}.")
            exit(1)
else:
    if timesteps < 1:
        print(f"Error: Given {timesteps} timesteps, require int > 1.")
        exit(1)
sim.run_timesteps(timesteps)