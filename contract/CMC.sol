// SPDX-License-Identifier: GPL-3.0
pragma solidity >=0.7.0 <0.9.0;
import "./ASC.sol";

contract CMC{
    ASC asc;
    struct credential{
        bytes32 Tx;
        bytes32 taux;
        bytes32 attrSets;
        address issuer;
        uint8 credentialType;
        bool Tlsb;
        bool taulsb;
        bool status;
    }
    address owner;
    mapping (bytes32 => credential) public credentials;
    event LogCredentialIssue(bytes32, address);
    event LogCredRevocation(bytes32);
    event LogCredRecover(bytes32);
    
    constructor(address a) {
        owner = a;
    }
    function setasc(address ascaddr) public {
        require (msg.sender == owner);
        asc = ASC(ascaddr);
    }


    function CredentialIssue(bytes32  Tx, bytes32 taux, bytes32 attrSets, uint8 CredType, bool  Tlsb, bool taulsb)public returns (bytes32) {
        //bytes32 tokenId = sha256(abi.encodePacked(_name,_value,msg.sender,_hash,now));
        assert(asc.setQuery(attrSets));
        bytes32 CredId = sha256(abi.encodePacked(Tx,taux,msg.sender,Tlsb, taulsb));
        credentials[CredId].Tx = Tx;
        credentials[CredId].taux = taux;
        credentials[CredId].issuer = msg.sender;
        credentials[CredId].attrSets = attrSets;
        credentials[CredId].credentialType = CredType;
        credentials[CredId].status = true;
        credentials[CredId].Tlsb = Tlsb;
        credentials[CredId].taulsb = taulsb;
        emit LogCredentialIssue(CredId, msg.sender);
        return CredId;
    }
    function CredRevocation(bytes32 CredId)public {
        if (msg.sender == credentials[CredId].issuer){
            credentials[CredId].status = false;
            emit LogCredRevocation(CredId);
        }
    }
    function CredRecover(bytes32 CredId)public {
        if (msg.sender == credentials[CredId].issuer){
            credentials[CredId].status = true;
            emit LogCredRecover(CredId);
        }
    }

}