set +e
mkdir -p /tutorial/
[ -f /home/snn.net ] && mv /home/snn.net /tutorial/
[ -f /home/arch.yaml ] && mv /home/arch.yaml /tutorial/
[ -f /home/loihi.yaml ] && mv /home/loihi.yaml /tutorial/
[ -f /home/dvs_challenge.npz ] && mv /home/dvs_challenge.npz /tutorial/
[ -f /home/dvs_challenge.py ] && mv /home/dvs_challenge.py /tutorial/
[ -f /home/TUTORIAL.md ] && mv /home/TUTORIAL.md /tutorial/
echo "SANA-FE Docker Image started successfully"